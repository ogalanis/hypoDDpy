#!/usr/bin/env python
# -*- coding: utf-8 -*-"
"""
Class handling the compilation of HypoDD.

    * "working_dir"/bin/hypoDD
    * "working_dir"/bin/ph2dt
    * "working_dir"/bin/hypoDD.inc

hypoDD and ph2dt are the binaries for the respective programs and hypoDD.inc is
the hypoDD.inc file used for the compilation.

If all three files are present and the hypoDD.inc that would be used for a new
compilation is identical to the one already present nothing will happen as the
end result would be the same.
"""
# import md5
import hashlib
import os
import shutil
import subprocess
import tarfile


# Specify the HypoDD version to be compiled.
HYPODD_ARCHIVE = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "src", "HYPODD_2.1b.tar.gz")
)

# Note this Hash is from the tar.gz file you got, either change this to your
# correct one using the function call md5.md5(open_file.read()).hexdigest()
# or simply comment out the line: if md5_hash != HYPODD_MD5_HASH:
# WCC added the hash pour his version and changed test from "!=" to "not in"
HYPODD_MD5_HASHES = [
    "ac7fb5829abef23aa91f1f8a115e2b45",
    "94228305b2370c4f3371fc6cb76f92c5",
]


class HypoDDCompilationError(Exception):
    """
    Exception that will be raised if anything during the compilation does not
    occur as planned.
    """

    pass


class HypoDDCompiler(object):
    """
    Class handling the HypoDD compilation.

    Usage
    =====

    >>> hyp_comp = HypoDDCompiler("temp_dir")
    >>> hyp_comp.configure()
    >>> hyp_comp.make()
    """

    def __init__(self, working_dir, log_function):
        """
        :param working_dir: The working directory. Everything will happen in
            there.
        :param log_function: Function to use to log activity.
        """
        # Set the log function.
        self.log = log_function
        # Set the working dir and create it if necessary.
        self.working_dir = working_dir
        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)
        # Make sure the given HypoDD archive is valid.
        self.verify_archive()
        # Setup and determine all the necessary paths.
        self.determine_paths()
        self.is_configured = False

    def verify_archive(self):
        """
        Method that checks if the HypoDD archive exists and that its md5 has is
        valid.
        """
        if not os.path.exists(HYPODD_ARCHIVE):
            msg = "HypoDD archive file could not be found"
            raise HypoDDCompilationError(msg)
        # Check if the file is correct.
        with open(HYPODD_ARCHIVE, "rb") as open_file:
            md5_hash = hashlib.md5(open_file.read()).hexdigest()
        # if md5_hash not in HYPODD_MD5_HASHES:
        #     msg = "md5 hash of the HypoDD archive is not correct"
        #     raise HypoDDCompilationError(msg)

    def determine_paths(self):
        self.paths = {}
        # Binary dir.
        self.paths["binary_dir"] = os.path.join(self.working_dir, "bin")
        if not os.path.exists(self.paths["binary_dir"]):
            os.makedirs(self.paths["binary_dir"])
        # Output files.
        self.paths["hypoDD_binary"] = os.path.join(
            self.paths["binary_dir"], "hypoDD"
        )
        self.paths["ph2dt_binary"] = os.path.join(
            self.paths["binary_dir"], "ph2dt"
        )
        # The hypoDD.inc files produced by any potential previous runs. After
        # the run, the currently used hypoDD.inc file will be copied there.
        self.paths["old hypoDD.inc file"] = os.path.join(
            self.paths["binary_dir"], "hypoDD.inc"
        )
        # Where to unpack the archive.
        self.paths["hypodd_unpack_dir"] = os.path.join(
            self.working_dir, "hypodd_src"
        )
        # Some paths in the unpacked archive.
        self.paths["make_directory"] = os.path.join(
            self.paths["hypodd_unpack_dir"], "HYPODD", "src"
        )
        # The resulting binaries directly after the compilation.
        self.paths["compiled_hypodd_binary"] = os.path.join(
            self.paths["make_directory"], "hypoDD", "hypoDD"
        )
        self.paths["compiled_ph2dt_binary"] = os.path.join(
            self.paths["make_directory"], "ph2dt", "ph2dt"
        )
        # The include directory.
        self.paths["include_dir"] = os.path.join(
            self.paths["hypodd_unpack_dir"], "HYPODD", "include"
        )
        # The hypoDD.inc file
        self.paths["hypoDD.inc"] = os.path.join(
            self.paths["include_dir"], "hypoDD.inc"
        )

    def configure(
        self,
        MAXEVE=3000,
        MAXDATA=2800000,
        MAXEVE0=50,
        MAXDATA0=60000,
        MAXLAY=30,
        MAXSTA=2000,
        MAXCL=200,
    ):
        """
        Configure the compilation.

        **hypoDD.inc configuration**

        The following parameters are used to configure the hypoDD.inc file. The
        default values are suitable for a medium sized problem. Adjust them if
        necessary.

        :param MAXEVE: Max number of events (must be at least the size of the
            number of events listed in the event file)
            Defaults to 3000.
        :param MAXDATA: Max number of observations (must be at least the size
            of the number of observations).
            Defaults to 2800000.
        :param MAXEVE0: Max number of events used for SVD. If only LSQR is
            used, MAXEVE0 can be set to 2 to free up memory.
            Defaults to 50.
        :param MAXDATA0: Max number of observations used for SVD. If only LSQR
            is used, MAXDATA0 can be set to 1 to free up memory.
            Defaults to 60000.
        :param MAXLAY: Max number of model layers.
            Defaults to 30.
        :param MAXSTA: Max number of stations.
            Defaults to 2000.
        :param MAXCL: Max number of clusters allowed.
            Defaults to 200.
        """
        # Set the hypodd_inc configuration.
        self.hypodd_inc_config = {
            "MAXEVE": MAXEVE,
            "MAXDATA": MAXDATA,
            "MAXEVE0": MAXEVE0,
            "MAXDATA0": MAXDATA0,
            "MAXLAY": MAXLAY,
            "MAXSTA": MAXSTA,
            "MAXCL": MAXCL,
        }

        self.is_configured = True

    def unpack_archive(self):
        """
        Unpacks the HypoDD archive to the hypodd_src subfolder in the working
        directory.
        """
        self.log("Unpacking HypoDD archive ...")

        unpack_dir = self.paths["hypodd_unpack_dir"]
        if os.path.exists(unpack_dir):
            shutil.rmtree(unpack_dir)
        os.makedirs(unpack_dir)

        tar = tarfile.open(HYPODD_ARCHIVE, "r:gz")
        tar.extractall(unpack_dir)
        self.log("Unpacking HypoDD archive done.")

    def make(self):
        if self.is_configured is not True:
            msg = "Compiler object need to be configured first."
            raise HypoDDCompilationError(msg)
        # Unpack the archive.
        self.unpack_archive()
        # Create the hypoDD_inc file.
        self.hypodd_inc_file = self.create_hypoDD_inc_file()
        # Check the current HypoDD compilation (if any).
        if self.is_current_hypodd_compilation_valid() is True:
            shutil.rmtree(self.paths["hypodd_unpack_dir"])
            self.log("Current compilation is up to date.")
            return
        # Finally compile it.
        self.compile_hypodd()
        # Cleanup.
        shutil.rmtree(self.paths["hypodd_unpack_dir"])

    def create_hypoDD_inc_file(self):
        """
        HypoDD uses static allocation and thus oftentimes has to be recompiled
        to suit a new problem size. The hypoDD.inc file is usually the only
        file that has to be changed. This is handled in this method.

        The default parameters are suitable for a medium sized problem.

        :return: A string containing the whole file.

        Original documentation:
        hypoDD.inc: Stores parameters that define array dimensions in hypoDD.
            Modify to fit size of problem and available computer memory.  If 3D
            raytracing is used, also set model parameters in vel3d.inc.

        Parameter Description:
        MAXEVE:   Max number of events (must be at least the size of the number
                  of events listed in the event file)
        MAXDATA:  Max number of observations (must be at least the size of the
                  number of observations).
        MAXEVE0:  Max number of events used for SVD. If only LSQR is used,
                  MAXEVE0 can be set to 2 to free up memory.
        MAXDATA0: Max number of observations used for SVD. If only LSQR is
            used, MAXDATA0 can be set to 1 to free up memory.
        MAXLAY:   Max number of model layers.
        MAXSTA:   Max number of stations.
        MAXCL:    Max number of clusters allowed.
        """
        # Do not mess with the indentation as it is important for Fortran77.
        hypoDD_inc = """
      integer*4 MAXEVE, MAXLAY, MAXDATA, MAXSTA, MAXEVE0, MAXDATA0
      integer*4 MAXCL
      parameter(MAXEVE   = {MAXEVE},
     &          MAXDATA  = {MAXDATA},
     &          MAXEVE0  = {MAXEVE0},
     &          MAXDATA0 = {MAXDATA0},
     &          MAXLAY   = {MAXLAY},
     &          MAXSTA   = {MAXSTA},
     &          MAXCL    = {MAXCL})""".format(
            MAXEVE=self.hypodd_inc_config["MAXEVE"],
            MAXDATA=self.hypodd_inc_config["MAXDATA"],
            MAXEVE0=self.hypodd_inc_config["MAXEVE0"],
            MAXDATA0=self.hypodd_inc_config["MAXDATA0"],
            MAXLAY=self.hypodd_inc_config["MAXLAY"],
            MAXSTA=self.hypodd_inc_config["MAXSTA"],
            MAXCL=self.hypodd_inc_config["MAXCL"],
        )
        # Remove the leading empty line.
        hypoDD_inc = hypoDD_inc[1:]
        return hypoDD_inc

    def is_current_hypodd_compilation_valid(self):
        """
        Returns True if the current compilation is ok, False otherwise. False
        should always trigger a new compilation.
        """
        # If the binary dir does not exist return False.
        if not os.path.exists(self.paths["binary_dir"]):
            return False
        # If the three file do not exist return False.
        if (
            not os.path.exists(self.paths["hypoDD_binary"])
            or not os.path.exists(self.paths["ph2dt_binary"])
            or not os.path.exists(self.paths["old hypoDD.inc file"])
        ):
            return False
        # Check if the newly created hypoDD.inc file is identical to the old
        # one.
        with open(self.paths["old hypoDD.inc file"], "r") as open_file:
            old_hypodd_file = open_file.read()
        if old_hypodd_file != self.hypodd_inc_file:
            return False
        return True

    def compile_hypodd(self):
        """
        Actually compiles HypoDD.
        """
        # Replace hypoDD.inc file with the custom one.
        os.remove(self.paths["hypoDD.inc"])
        with open(self.paths["hypoDD.inc"], "w") as open_file:
            open_file.write(self.hypodd_inc_file)
        # Compile it.
        self.log("Compiling HypoDD ...")
        sub = subprocess.Popen(
            "make",
            cwd=self.paths["make_directory"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
        )
        self.log(sub.stdout.read())
        retcode = sub.wait()
        if retcode != 0:
            msg = "Problem compiling HypoDD."
            raise HypoDDCompilationError(msg)
        # Check if the output files have been created.
        if not os.path.exists(
            self.paths["compiled_hypodd_binary"]
        ) or not os.path.exists(self.paths["compiled_ph2dt_binary"]):
            msg = "The binary output files could not be found."
            raise HypoDDCompilationError(msg)
        # Move the binary files and the hypoDD.inc file.
        shutil.move(
            self.paths["compiled_hypodd_binary"], self.paths["hypoDD_binary"]
        )
        shutil.move(
            self.paths["compiled_ph2dt_binary"], self.paths["ph2dt_binary"]
        )
        shutil.move(
            self.paths["hypoDD.inc"], self.paths["old hypoDD.inc file"]
        )
        self.log("Compiling HypoDD done.")
