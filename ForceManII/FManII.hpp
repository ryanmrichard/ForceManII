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
/** \file FManII.hpp
 * 
 * \brief This header is meant to be the main API to ForceManII.  Anything users
 *     need to interface with ForceManII goes here.
 * 
 * \version 0.1
 * \date October 14, 2016 at 10:17 PM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FMANII_FMANII_HPP
#define FMANII_FMANII_HPP

#include "ForceManII/FManIIDefs.hpp"
#include "ForceManII/ForceField.hpp"
#include "ForceManII/InternalCoordinates.hpp"
#include <istream>
#include <cmath>

///Namespace for all code associated with ForceManII
namespace FManII {

///An array of internal coordinates arranged by type
using CoordArray=std::map<IntCoord_t,std::unique_ptr<IntCoords>>;


/**\brief Given a force field file in Tinker format makes a ForceField object
 *
 * Optionally one may specify their own unit conversions as well
 * \param[in] file an istream instance loaded with a Tinker formated string
 * \param[in] kcalmol2au The conversion from kcal/mol to Hartrees
 * \param[in] ang2au The conversion from Angstroms to Bohr
 * \param[in] deg2rad The conversion from degrees to radians
 * \return Your parsed force field
 */
ForceField parse_file(std::istream&& file,
                      double kcalmol2au=1.0/627.5096,
                      double ang2au=1.889725989,
                      double deg2rad=M_PI/180.0);


/**\brief A function that processes the input and returns a set of objects set
 *     up for use in FManII
 * 
 * See the documentation for each type for information on how to set it up.
 * Each type is actually just a typedef of STL containers so no other resources
 * of FManII are needed to run this function.
 * 
 * \param[in] Carts The Cartesian coordinates in a.u. of each atom in the form
 *                  \f$x,y,z\f$ of atom 1, then \f$x,y,z\f$ of atom 2, etc.
 * \param[in] Types An array of the atom types of each atom
 * \param[in] Params The parameters for the force field
 * \param[in] Conns The connectivity information
 * \param[in] chg14scale How much should 1-4 electrostatics be scaled by
 * \param[in] vdw14scale How much should 1-4 VDW terms be scaled by
 * 
 */
CoordArray get_coords(const std::vector<double>& Carts,
                     const AtomTypes& Types,
                     const ParamTypes& Params,
                     const ConnData& Conns,
                     double chg14scale=1.0,double vdw14scale=1.0);
} //End namespace FManII

#endif /* End header guard */

