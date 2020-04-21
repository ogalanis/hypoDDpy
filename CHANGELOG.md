hypoddpy/__init__.py
    Changed "from hypodd_relocator import HypoDDRelocator"
    to "from .hypodd_relocator import HypoDDRelocator"
ran 2to3 -p -v -w
changed "import md5" to "import hashlib"
changed hypodd_relocator.HyopDDRelocator._parse_station_files() to
    read StationXML files
Set subrocess.Popen(unversal_newlines=True) in hypodd_compiler.compile_hypodd()
    so that the stdout output will be a text string (as in Python 2)
    rather than a byte string
modified station_id to not include network if {net}.{sta} > 7 characters
2 to 3 bug?: changed Exception.message to str(Exception) (lines 1104)
**MUST DO SOMETHING ABOUT NEGATIVE ELEVATIONS (HYPODD DOESNT HANDLE!)**
# Probably need to ship stations up so deepest is at zero, then
# shift model up by as much, then go backwards when reading in results
added shift_stations attribute (and associated code) to HypDDRelocator class
  * Still needs to shift input and output events

- corrected bug in  _write_ph2dt_inp_file(self) where maxsep was 
  calculated using depth differences in meters instead of km
- HAND-CHANGED iteration distance values (should be calculable or 
  enterable for DATA_WEIGHTING_AND_REWEIGHTING variable in
  _write_hypoDD_inp_file().  (were 0.1, 0.005; now 5, 2.5)
- hypodd_relocator.setup_velocity_model(): Added check for > 30 model layers
- hypodd_relocator.compile_hypodd(): Changed MAXDATA to 3000000 (should
  be configurable)
- hypodd_relocator._create_output_event_file(): changed res_id.getRefferedObject()
    to res_id.get_referred_object()