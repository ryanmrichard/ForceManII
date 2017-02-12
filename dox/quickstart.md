# Quick Start               {#quickstart}
========================

This page is intended to help you quickly add the ForceManII project to your
existing code.

## Download and Compile
The official repository for ForceManII is on  GitHub at:
[https://github.com/ryanmrichard/ForceManII]
(https://github.com/ryanmrichard/ForceManII).
First you need to download the source using the normal git commands:

```.sh
git clone https://github.com/ryanmrichard/ForceManII.git <where_to_put_it>
```
The last argument is optional, if you don't use it the source will be downloaded
into a subdirectory entitled `ForceManII` of the directory in which you ran git.
For this tutorial we assume you maintained the default name.  Once it is done
downloading:

```.sh
cd ForceManII
#Configure the build
cmake -H. -Bbuild -DCMAKE_CXX_COMILER=<path_to_your_c++_compiler> \
                  -DCMAKE_INSTALL_PREFIX=<where_to_put_the_final_library>

#Build the library
cd build && make

#Optionally test the library
ctest

#Install the library (may require sudo depending on install location)
make install

#Enjoy science
```
ForceManII has no dependencies aside from a C++11 compliant compiler and a
relatively modern version of CMake (circa 2014 or later).  That being said on
most systems you will not even have to specify the CXX compiler as CMake will
detect one automatically.  As a CMake project,
ForceManII strives to honor the usual CMake variables where appropriate;
therefore power users should feel free to pass additional CMake variables for
more fine-grained control.

## Using the ForceManII API

The absolute simplest call to ForceManII is:

```.cpp
#import <ForceManII/FManII.hpp>

auto deriv=
  FManII::run_forcemanii(order,carts,conns,FManII::get_ff(ff_name),types);
```
Here:

- `order` is the derivative order you want (0=energy,1=gradient,2=Hessian,
*etc.*).
- `carts` is a 3 by number of atoms `std::vector<double>` where
`carts[i*3+j]` is the \f$j\f$-th Cartesian component (\f$j\f$=0 is \f$x\f$,
\f$j\f$=1 is \f$y\f$, \f$j\f$=2 is \f$z\f$) of the \f$i\f$-th atom.
- `conns` is an `std::vector<std::vector<size_t>>` (basically a matrix where the
rows have variable lengths) such that `conns[i][j]` is the index of the
\f$j\f$-th atom connected to the \f$i\f$-th atom (order being the same as carts)
- FManII::get_ff() is a helper function that will return any of the built-in
force fields given...
- `ff_name` the name of the force field as an `std::string` currently we support
  - `AMBER99`
  - `CHARMM22`
  - `OPLSAA`
- `types` is an `std::vector<size_t>` where `types[i]` is the atom type of atom
\f$i\f$.  Numbering of atom types varies from force field to force field look in
the ForceFields directory for more info
- `deriv` will be an `std::map<std::string,std::vector<double>>` where the key
  is a descriptive name of what derivative component you are looking at, *e.g.*
  the harmonic bond stretching term, and the value is the derivative in C++
  order, *i.e.* row-major

Where applicable, all units are atomic units, *i.e.* derivatives are Hartrees
over Bohrs to the derivative order, input Cartesian coordinates are in Bhors.

The API to ForceManII is designed to be as flexible as possible while still
maintaining simplicity.  Therefore much of the flow of the program can be
controlled from outside the library.  For example, say you wanted to use a force
field that is not included in ForceManII, you can do this by:

~~~.cpp
//Conversions are optional, default values exist
FManII::ForceField my_ff=FManII::parse_file
(
    std::move(std::istream("path/to/.prm/file")),
    your_kcalmol_2_hartree_conversion,
    your_angstrom_2_bohr_conversion,
    your_degree_2_radian_conversion
);
~~~

You would now replace FManII::get_ff()`with `my_ff` in the call to
FManII::run_forcemanii().

The FManII::ForceField object is relatively simple, so if you wanted to make
your own ForceField all you would need to do is set the membere appropriately
and use the resulting instance.

The FManII::run_forcemanii() function is actually a thin wrapper around a series
of steps: compute internal coordinates, assign parameters, and compute the
derivative.  The ForceManII API allows you to call each of those functions
manually so you can further customize the command, or bypass a command
completely and just use your own objects.


