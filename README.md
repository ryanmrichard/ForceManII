# ForceManII: An MM code for QM people

TODO: CI integration

User and Developer manual is available [here](https://ryanmrichard.github.io/ForceManII/).

## History
ForceManII was originally a code written to replace Q-Chem's original MM code, 
which is entitled ForceMan (short for force-field manager).  I (Ryan M. Richard
), started work on ForceManII as a graduate student at Ohio State where I was
part of the Q-Chem development team.  I then went on to do a postdoc at Georgia
 Institue of Technology where I became a Psi4 developer.  In an effort to still
 maintain ForceManII the decision was made to make it a stand alone library that
 could be called from either Q-Chem or Psi4.

## Features

TODO: Explain to the world why this is the greatest MM code ever


## Compilation

ForceManII has no dependencies other than a C++11 compliant compiler.  The
following should thus be sufficient to compile it:

~~~.sh
#In the top-level directory of the ForceManII source you got from GitHub
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=path/to/your/c++/compiler
cd build && make

#Optionally run the test suite
ctest

#Install (may require root depending on install location)
make install
~~~

Note that due to the amount of information that the tests test compilation of
the tests may be a bit slow (on a semi-modern computer, the entire build,
including tests, should be less than 5 minutes though).

There are a couple of other useful CMake variables you may want to consider
setting:

- `CMAKE_INSTALL_PREFIX` : Sets the install location
- `CMAKE_BUILD_TYPE` : Defaults to `release`, but set to `debug` if you want
  debug symbols

As usual these would be appended to your cmake command with the `-D` flag,
*e.g.*, to change the install prefix your new CMake command would be:

~~~.sh
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=path/to/your/c++/compiler\
                  -DCMAKE_INSTALL_PREFIX=where/you/want/ForceManII/installed
~~~


## Quick-Start

If all you want to do is compute an MM energy for an existing force field this
quick-start guide should give you all you need.  Assuming you have ForceManII
compiled and ready to be used:

In the *.cpp* file you want to call ForceManII from:
~~~.cpp
#include<ForceManII/FManII.hpp>
//Optional, if you want to use one of the force fields that ships with FManII
//e.g. Amber99, include its header file
//#include <ForceFields/amber99.hpp>

double my_fxn(){
    //Conversions are optional, default values exist, skip this if you are
    //using a force field that ships with FManII
    FManII::ForceField my_ff=FManII::parse_file
    (
        std::move(std::istream("path/to/.prm/file")),
        your_kcalmol_2_hartree_conversion,
        your_angstrom_2_bohr_conversion,
        your_degree_2_radian_conversion
    );

    //You somehow need to get the coordinates of your system, in bohr, into
    //an NAtoms by 3 long vector
    std::vector<double> my_carts=magic_function_that_gets_coordinates();

    //You somehow need to get the atom types of the atoms in your system
    //AtomTypes is a typedef of std::vector<size_t>
    FManII::AtomTypes atom_types=magic_function_that_gets_types();

    //You somehow need to get the connectivity of your system
    //This is an NAtoms long list where element i is a list of the atoms
    //connected to atom i

    //ConnData is a typedef of std::vector<std::vector<size_t>>
    FManII::ConnData conns=magic_function_that_gets_conns();

    //Get the internal coordinates of your system
    FManII::CoordArray Coords=FManII::get_coords
    (
        my_carts,
        atom_types,
        conns,
        my_ff //would be FManII::amber99 if you were using included Amber99 ff
     );

    //Computes the requested derivative
    return FManII::compute_deriv(CoordArray,derivative_order);

}
~~~

As you can hopefully see the API design is designed to allow the user to
customize each step of the process.

- The decoupling of the three main steps allows the user to modify input to the
  next step without modifying ForceManII
- `parse_file` is designed to obtain the details of a force field
   - Force field files should be in Tinker format.
   - Files can be custom force fields or modified/extended versions of the
     force fields that ship with ForceManII
   - One can also fill a `FManII::ForceField` instance manually for a more
     custom experience
- `get_coords` determines the internal coordinates of a system and assigns
   params
   - Internal coordinates *e.g.* bonds, angles, torsions, etc.
   - Combined determine/assign procedure is for effeciency
- `compute_deriv` computes the requested derivative
  - note that the 0-th derivative is the energy, 1-st is the gradient, *etc.*
