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

#ifndef FMANII_FMANII_HPP
#define FMANII_FMANII_HPP

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

///Convenience classes for common model/coord choices
///@{
struct HarmonicBond:public FFTerm{HarmonicBond();};

struct HarmonicAngle:public FFTerm{HarmonicAngle();};

struct FourierTorsion:public FFTerm{FourierTorsion();};

struct FourierImproperTorsion:public FFTerm{FourierImproperTorsion();};

struct LJ14:public FFTerm{LJ14();};

struct LJPair:public FFTerm{LJPair();};

struct Electrostatics14:public FFTerm{Electrostatics14();};

struct ElectrostaticsPair:public FFTerm{ElectrostaticsPair();};

///@}


///Available hard-coded force fields
extern const ForceField amber99;

/**\brief Given a force field file in Tinker format makes a ForceField object
 *
 * Optionally one may specify their own unit conversions as well
 * \param[in] file an istream instance loaded with a Tinker formated string
 * \param[in] kcalmol2au The conversion from kcal/mol to Hartrees
 * \param[in] ang2au The conversion from Angstroms to Bohr
 * \param[in] deg2rad The conversion from degrees to radians
 * \return Your parsed force field in atomic units
 */
ForceField parse_file(std::istream&& file,
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
CoordArray get_coords(const Vector& Carts,
                      const ConnData& Conns);


/**\brief A function that assigns the final parameters to a system
 *
 *
 */
ParamSet assign_params(const CoordArray& coords,
                       const ForceField& ff,
                       const IVector& types);

DerivType deriv(size_t order,
                const ForceField& ff,
                const ParamSet& ps,
                const CoordArray& coords);

inline DerivType run_forcemanii(size_t order,
                                const Vector& Carts,
                                const ConnData& conns,
                                const ForceField& ff,
                                const IVector& types){
    const CoordArray coords=get_coords(Carts,conns);
    return deriv(order,ff,assign_params(coords,ff,types),coords);

}


} //End namespace FManII

#endif /* End header guard */

