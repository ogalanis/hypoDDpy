#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
hypoddpy - Simple HypoDD workflow controller.

hypoddpy is able to do a full relocation using HypoDD (Waldhauser, 2000). The
only inputs required are the event information as QuakeML, the station
information in SEED/XSEED and the waveform files in any format ObsPy can read.

The output is a QuakeML file that is identical to the input event information
and additionally every event that has been relocated will have an addtional
Origin.


[Waldhauser, F., & Ellsworth, W. L. (2000). A Double-Difference Earthquake
Location Algorithm: Method and Application to the Northern Hayward Fault,
California. Bulletin of the Seismological Society of America, 90(6), 1353-1368.
doi:10.1785/0120000006]

:copyright:
    lion krischer (krischer@geophysik.uni-muenchen.de), 2012
:license:
    gnu lesser general public license, version 3
    (http://www.gnu.org/copyleft/lesser.html)
"""
from __future__ import with_statement

from setuptools import setup
import os

LOCAL_PATH = os.path.abspath(os.path.dirname(__file__))
DOCSTRING = __doc__.split("\n")

NAME = 'hypoddpy'
AUTHOR = 'Lion Krischer'
AUTHOR_EMAIL = 'krischer@geophysik.uni-muenchen.de'
URL = 'None so far...'
LICENSE = 'GNU General Public License, version 3 (GPLv3)'
KEYWORDS = ['seismology', 'earthquakes', 'relocation']
INSTALL_REQUIRES = ['obspy', 'progressbar']
ENTRY_POINTS = {}


def getVersion():
    """
    Get the current version of the module.
    """
    version_file = os.path.join(LOCAL_PATH, 'hypoddpy', 'VERSION.txt')
    with open(version_file) as f:
        version = f.read().strip()
    return version


def setupPackage():
    # setup package
    setup(
        name=NAME,
        version=getVersion(),
        description=DOCSTRING[1],
        long_description="\n".join(DOCSTRING[3:]),
        url=URL,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        license=LICENSE,
        platforms='OS Independent',
        keywords=KEYWORDS,
        packages=['hypoddpy'],
        package_dir={'hypoddpy': 'hypoddpy'},
        zip_safe=False,
        install_requires=INSTALL_REQUIRES,
        entry_points=ENTRY_POINTS,
    )


if __name__ == '__main__':
    setupPackage()
