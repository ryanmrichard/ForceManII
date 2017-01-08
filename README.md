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
- Reliance on atomic units
- A library like model
  - Most MM packages are just that packages.  They are not designed to be called
    from other programs.  This makes it challenging to write QM/MM models or
    use MM results in traditional QM packages.
- Easy addition of new force fields with (possibly) new functional forms
  - This requires a flexible API, which is described in the documentation
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

TODO: Make a CMake config file for ForceManII to facilitate adding it to other
CMake projects

## Quick-Start

If all you want to do is compute an MM energy for an existing force field this
quick-start guide should give you all you need.  Assuming you have ForceManII
compiled and ready to be used:

In the *.cpp* file you want to call ForceManII from:
~~~.cpp
#include<ForceManII/FManII.hpp>

void my_fxn(){
    //Step 1: Get some information that you should have laying around such as...

    //...the Cartesian coordinates of your system, in bohr, such that
    //my_carts[i*3+j] is the j-th Cartesian coordinate of atom i
    std::vector<double> my_carts=magic_function_that_gets_coordinates();

   //...a vector of each atom's atom type such that atom_types[i] is the
   //atom type of atom i
   std::vector<size_t> atom_types=magic_function_that_gets_types();

   //...and the connectivity of your system, such that conns[i] is an
   //std::set<size_t> of the atoms connected to atom i
   std::vector<std::set<size_t>> conns=magic_function_that_gets_conns();

   //Step 2: Obtain a ForceField instance either by...

   //...parsing an existing Tinker-style parameter file
   FManII::ForceField my_ff=FManII::parse_file
    (
        std::move(std::istream("path/to/.prm/file")),
        your_kcalmol_2_hartree_conversion,
        your_angstrom_2_bohr_conversion,
        your_degree_2_radian_conversion
    );

    //...using one of the included, hard-coded force fields
    //my_ff=FManII::amber99 //Uncomment this line for AMBER99

    //...or filling one in yourself
    //my_ff=magic_fucntion_that_makes_a_custom_force_field();

    //Step 3: Get your requested derivative
    FManII::DerivType deriv=
        FManII::run_fmanii(0,my_carts,conns,my_ff,atom_types);

    //Step 4: Do science!!!!
}
~~~

A couple of notes on the above to hopefully fill in the missing details:
- All atoms must have a type
  - It is our hope to automate this step based on the connectivity in an
    upcoming release
- A future release will be able of automatically computing the connectivity
- `FManII::run_fmanii` is actually a wrapper around several functions:
  - A function that finds the internal coordinates
  - A function that assigns parameters
  - And a function that computes the derivative using those coordinates and
    derivatives
  - The user can skip any of these functions and supply their own information
    for an even more customizable experience
- The type of `FManII::DerivType` is
  `std::map<std::pair<std::string,std::string>,std::vector<double>>`.  I know
  this looks scary, but it's really quite simple if we break it down.
  - The key of the map is of type `std::pair<std::string,std::string>`,which is
    a type we call FManII::FFTerm.  It is just a pair of names describing the
    model (harmonic oscillator, Fourier series, etc.) and the internal
    coordinate type (bond, angle, torsion, etc.)
    - *e.g.* a harmonic-angle bending term is
      `FFTerm({"HARMONICOSCILLATOR","ANGLE"})` and a harmonic-bond stretching
       term is `FFTerm({"HARMONICOSCILLATOR","BOND"})`
  - The value is the actual derivative as a vector (so the energy is a single
    element vector, the gradient is a 3 by natoms long vector, and the Hessian
    is a 3 times natoms by 3 times natoms matrix flattened into a vector)
  - For example to print the energies try:
    ~~~.cpp
    for(auto name_deriv:deriv){
        std::string &model_name=name_deriv.first.first;
        std::string &coord_name=name_deriv.first.second;
        std::vector<double> &derivi=name_deriv.second;
        std::cout<<"The energy of the "<<coord_name<<"s is "<<derivi[0]<<
           " (a.u.) if we assume they are modeled with a "<<model_name<<
           " potential."<<std::endl;
    }
    ~~~
