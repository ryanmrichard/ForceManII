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
/** \file FManIIDefs.hpp
 *
 * \brief This header contains enums and other defs used
 * throughout the code
 *
 * \version 0.1
 * \date December 10, 2016 at 2:38 PM (EST)
 *
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 *
 * Additional contributions by:
 *
 */
#ifndef FMANIIDEFS_HPP
#define FMANIIDEFS_HPP

#include <map>
#include <vector>
#include <array>
#include <memory>

namespace FManII {

///These are the recognized types of parameters
enum Param_t {
    K,///<A force constant
    r0,///<The equilibrium value
    amp,///<The amplitude for Fourier series
    phi,///<The phase shift for Fourier series
    n,///<The periodicity for Fourier series
    amp2,///<The amplitude for a two-part Fourier series
    phi2,///<The phase shift for a two-part Fourier series
    n2,///<The periodicity for a two-part Fourier series
    amp3,///<The amplitude for a three-part Fourier series
    phi3,///<The phase shift for a three-part Fourier series
    n3,///<The periodicity for a three-part Fourier series
    q,///<The charge, in a.u., for point-charge, point-charge
    sigma,///<The minimum diameter of a 6-12 potential
    epsilon,///<The well depth of a 6-12 potential
};

///These are the recognized types of IntCoords
enum IntCoord_t {
    BOND,///<A bond
    UBPAIR,///<A 1,3 pair
    PAIR,///<A pair that is not a bond or a 1,3 pair
    ANGLE,///<An angle
    TORSION,///<A torsion
    IMPTORSION,///<An improper torsion angle
    ELECTROSTATICS,///<A charge-charge interaction
    LENNARD_JONES,///<A 6-12 Lennard-Jones potential
};

///These are the recognized combination rules
enum CombRule_t {
    ARITHMETIC,///<A normal average \f$\frac{1}{N}\sum_{i=1}^Nx_i\f$
    GEOMETRIC,///<Geometric average \f$\left(\prod_{i=1}^Nx_i\right)^{1/N}\f$
};

///Array where element i respectively is the class type and atom type of atom i
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

}


#endif // FMANIIDEFS_HPP
