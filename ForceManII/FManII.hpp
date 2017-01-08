/*
 * Copyright (C) 2016 Ryan M. Richard <ryanmrichard1 at gmail.com>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */

#pragma once

#
#include "ForceManII/FManIIDefs.hpp"
#include "ForceManII/ForceField.hpp"
#include "ForceManII/InternalCoordinates.hpp"
#include "ForceManII/ModelPotential.hpp"
#include "ForceManII/FFTerm.hpp"
#include "ForceManII/ModelPotentials/HarmonicOscillator.hpp"
#include "ForceManII/ModelPotentials/LennardJones.hpp"
#include "ForceManII/ModelPotentials/FourierSeries.hpp"
#include "ForceManII/ModelPotentials/Electrostatics.hpp"

#include <istream>
#include <cmath>

///Namespace for all code associated with ForceManII
namespace FManII {

///Available hard-coded force fields
extern const ForceField amber99;
extern const ForceField oplsaa;
extern const ForceField charmm22;

///Convenience functions for making model potentials, internal coordinates, and
///force field terms
///@{

///Returns the model potential associated with the given key
std::shared_ptr<ModelPotential> get_potential(const std::string& name);

///Returns the internal coordinate associated with the given key
std::shared_ptr<InternalCoordinates> get_intcoord(const std::string& name);

///Makes a force field term with given model and internal coordinate
inline FFTerm get_term(const FFTerm_t& name){
    return FFTerm(get_potential(name.first),get_intcoord(name.second));
}
///@}

/**\brief Given a force field file in Tinker format makes a ForceField object
 *
 * Optionally one may specify their own unit conversions as well
 * \param[in] file an istream instance loaded with a Tinker formated string
 * \param[in] is_charmm Should improper torsions be interpreted as being in
 *              Tinker or CHARMM order
 * \param[in] kcalmol2au The conversion from kcal/mol to Hartrees
 * \param[in] ang2au The conversion from Angstroms to Bohr
 * \param[in] deg2rad The conversion from degrees to radians
 * \return Your parsed force field in atomic units
 */
ForceField parse_file(std::istream&& file,
                      bool is_charmm=false,
                      double kcalmol2au=1.0/627.5096,
                      double ang2au=1.889725989,
                      double deg2rad=M_PI/180.0);


/**\brief A function that processes the input and returns a set of internal
 *        coordinates
 *
 * \note Atom ordering refers to whatever order the user gives us the atoms in.
 *       This order is assumed arbitrary, but consistent (i.e. element i of
 *       \p Carts corresponds to the same atom as element i of \p Types, which
 *       corresponds to the same atom as element i of \p Conns)
 *
 *
 * \param[in] Carts The Cartesian coordinates in a.u. of each atom in the form
 *                  \f$x,y,z\f$ of atom 1, then \f$x,y,z\f$ of atom 2, etc.  If
 *                  you like this is an natoms by 3 matrix with atoms on the
 *                  rows and Cartesian coordinates on the columns.
 * \param[in] Types An array so that element i is the atom type of atom i
 * \param[in] Conns The connectivity information for your system such that
 *                  element i is a vector of atoms bonded to atom i (all atoms
 *                  must have a vector associated with them, even if it is
 *                  empty
 * \param[in] ff A ForceField instance with the details of your forcefield
 * \return Your system's internal coordinates, in a.u.
 * 
 */
Molecule get_coords(const Vector& Carts,
                      const ConnData& Conns);


/**\brief A function that assigns the final parameters to a system
 *
 *
 *  \param[in] coords The internal coordinates of the system
 *  \param[in] ff The force field to use for assigning parameters
 *  \param[in] types A map from atom number to atom type
 *  \param[in] skip_missing If true missing parameters will be counted as zero
 *  \return The parameters of your system
 *
 */
ParamSet assign_params(const Molecule& coords,
                       const ForceField& ff,
                       const IVector& types,
                       bool skip_missing=true);

DerivType deriv(size_t order,
                const ForceField& ff,
                const ParamSet& ps,
                const Molecule& coords);

inline DerivType run_forcemanii(size_t order,
                                const Vector& Carts,
                                const ConnData& conns,
                                const ForceField& ff,
                                const IVector& types){
    const Molecule coords=get_coords(Carts,conns);
    return deriv(order,ff,assign_params(coords,ff,types),coords);

}


} //End namespace FManII

