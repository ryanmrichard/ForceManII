# Quick Start               {#quickstart}
========================

The call sequence to the ForceManII library will be very similar in almost all
cases with some exceptions when the user desires more customability.

In general the call sequence will be:

1. `parse_file` to generate a ForceField instance
2. `get_coords` to generate a CoordArray instance for the system
3. `get_params` to parameterize each coordinate of the system
4. `deriv` to compute the requested derivative

For convenience we define a function `run_forcemanii` that does this sequence
for you returning the requested derivative.

Consequentially a typical usage of ForceManII is something like:

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


