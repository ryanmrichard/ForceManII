# ForceManII: An MM code for QM people

TODO: CI integration

User and Developer manual is available
[here](https://ryanmrichard.github.io/ForceManII/) and is the primary source of
information regarding usage of ForceManII (aside from this readme).

## Features

ForceManII is designed to facilitate running MM computations from QM codes. This
amounts to several important features in our opinion:

- Usage of data typically available in a QM code such as Cartesian coordinates
  and atomic numbers
  - At the moment we also need the connectivity of the atoms and the force field
  parameters for each atom, but we hope to remove these limitations in future
  releases
- A library like model
  - Most MM packages are just that packages.  They are not designed to be called
    from other programs.  This makes it challenging to write QM/MM models or
    use MM results in traditional QM packages.
- Easy addition of new force fields with (possibly) new functional forms
- Object-oriented design written in C++
- Minimal (currently no) dependencies


## Compilation

ForceManII has no dependencies other than a C++11 compliant compiler.  The
following should thus be sufficient to compile it:

~~~.sh
#In the top-level directory of the ForceManII source you got from GitHub
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=path/to/your/c++/compiler
cd build && make

#Optionally run the test suite (still in the build directory)
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

void my_fxn(){
    //Conversions are optional, default values exist, skip this call if you are
    //using a force field that ships with FManII and replace all instances of
    //my_ff with FManII::amber99
    FManII::ForceField my_ff=FManII::parse_file
    (
        std::move(std::istream("path/to/.prm/file")),
        your_kcalmol_2_hartree_conversion,
        your_angstrom_2_bohr_conversion,
        your_degree_2_radian_conversion
    );

    //You somehow need to get the coordinates of your system, in bohr, into
    //an NAtoms by 3 long vector
    //thus my_carts[i*3+j] is the j-th Cartesian coordinate of atom i
    std::vector<double> my_carts=magic_function_that_gets_coordinates();

    //You somehow need to get the atom types of the atoms in your system
    //AtomTypes is a typedef of std::vector<size_t>
    //thus atom_types[i] is the atom type of atom i
    FManII::AtomTypes atom_types=magic_function_that_gets_types();

    //You somehow need to get the connectivity of your system
    //This is an NAtoms long list where element i is a list of the atoms
    //connected to atom i

    //ConnData is a typedef of std::vector<std::vector<size_t>>
    //thus conns[i][j] is the j-th atom connected to atom i
    FManII::ConnData conns=magic_function_that_gets_conns();

    //Get the derivative you want
    //DerivType is a typedef of std::map<FManII::FFTerm_t,std::vector<double>>
    //thus `deriv.at({HARMONICOSCILLATOR,BOND})[i]` is the i-th component of the
    //harmonic bond-stretching term of the force field
    FManII::DerivType deriv=
        FManII::run_fmanii(0,my_carts,conns,my_ff,atom_types);
}
~~~

As you can hopefully see the API design is designed to allow the user to
customize each step of the process.

- The decoupling of the main steps allows the user to modify input to the
  next step without modifying ForceManII
- `parse_file` is designed to obtain the details of a force field
   - Force field files should be in Tinker format.
   - Files can be custom force fields or modified/extended versions of the
     force fields that ship with ForceManII
   - One can also fill a `FManII::ForceField` instance manually for a more
     custom experience
- `run_fmanii` computes the requested derivative
  - actually has three sub calls get_coords, get_params, deriv that
    respectively compute the internal coordinates, assign parmaeters to them,
    and then ultimately compute the derivative
    - These functions may be called by the user for even more fine grained
      control over the API
  - note that the 0-th derivative is the energy, 1-st is the gradient, *etc.*
