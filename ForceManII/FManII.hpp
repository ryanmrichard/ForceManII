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

#include <map>
#include <vector>
#include <array>
#include <memory>

///Namespace for all code associated with ForceManII
namespace FManII {

class IntCoords;

///These are the recognized types of parameters
enum Param_t {
    K,///<A force constant
    r0,///<The equilibrium value
    amp,///<The amplitude for Fourier series
    phi,///<The phase shift for Fourier series
    n///<The periodicity for Fourier series
};

///These are the recognized types of IntCoords
enum IntCoord_t {
    Bond,///<A bond
    UBPair,///<A 1,3 pair
    Pair,///<A pair that is not a bond or a 1,3 pair
    Angle,///<An angle
    Torsion,///<A torsion
    ImpTorsion,///<An improper torsion angle
};

///An array of internal coordinates arranged by type
using CoordArray=std::map<IntCoord_t,std::unique_ptr<IntCoords>>;

///Array where element i respectively is the atom type and VDW type of atom i
using AtomTypes=std::vector<std::array<size_t,2>>;

///Array such that element i is a vector of the atoms bonded to atom i
using ConnData=std::vector<std::vector<size_t>>;

/**\brief An object to hold the parameters of a force field
 *
 *  If you like you can think of this as a rank 3 tensor.  The first rank tells
 *  us which of the recognized internal coordinates the parameters are for.  The
 *  second rank tells us which of the recognized types of parameters the set is
 *  and the third rank tells us the atom type.  Putting that all together, the
 *  bond force constant between say atoms of type 3  and 4 is given by:
 *  \code
 *  ParamTypes params;//Assume we already had this
 *  double k=params[FManII::Bond][FManII::K][{3,4}];
 *  \endcode
 */
using ParamTypes=std::map<IntCoord_t,std::map<Param_t,
        std::map<std::vector<size_t>,double>>>;

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
 * 
 */
CoordArray get_coords(const std::vector<double>& Carts,
                     const AtomTypes& Types,
                     const ParamTypes& Params,
                     const ConnData& Conns);
} //End namespace FManII

#endif /* End header guard */

