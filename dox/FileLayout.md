ForceManII File Layout                                            {#file_layout}
======================

The purpose of this document is describe how the ForceManII source code is laid
out.  After reading it you should have some semblance of where I want files to
go and what is where.

At the top level should be the following directories:
- bin : Contains scripts for generating source files and test data
- dox : The source for the documentation
- ForceFields : A collection of force fields in Tinker format
- ForceManII : The main source directory
- tests : A set of unit tests

also in the top-level directory is the main readme for the GitHub page, the
project's license, a file detailing how to contribute, and the main build file.

## Directories and Files in More Detail

### bin

Particularly for the tests, there is a need to dump the correct answer to a
file.  For example, if you open `tests/ubiquitin_intcoords.cpp` you will find
a 761,830 line file that is every bond, angle, torsion, improper torsion, 1-4
pair, and pair in the protein ubiquitin.  Such files are crucial for unit tests
as they provide points of comparision that are unchanging for intermediate steps
in large computations.  Obviously, we don't want to generate such files by hand.
For this reason I have written a series of Python scripts that will generate
various types of data for you.

\warning The Python scripts are hastily written and not inteneded for production
they are highly messy (lots of copy paste code) and non-optimized.

The included scripts are:
- FManII.py : This is a sort of header file with functionality and definitions
  common to the other scripts
- ParseFF.py : This script generates the included hard coded force fields found
  in `ForceManII/ForceFields`
- FindIntCoords.py : This script generates the various `tests/*_intcoords.?pp`
  files used in `tests/TestIntCoords.cpp`
- AssignParameters.py : This script generates the various `tests/*_params.?pp`
  files used in `tests/TestAssignParams.cpp`
- MakeTest.py : This script generates the various `tests/<molecule_name>.?pp`
  files used in energy derivative checks.

\todo describe the remaining folders and files

### dox

### ForceFields

### ForceManII

### tests
