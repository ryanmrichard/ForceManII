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
    q,///<The charge, in a.u., for point-charge, point-charge
    sigma,///<The minimum diameter of a 6-12 potential
    epsilon,///<The well depth of a 6-12 potential
};

///These are the recognized types of IntCoords
enum IntCoord_t {
    BOND,///<A bond
    PAIR13,///<A 1,3 pair
    PAIR14,///<A 1,4 pair
    PAIR,///<A pair that is not a 1,2; 1,3; or 1,4 pair
    ANGLE,///<An angle
    TORSION,///<A torsion
    IMPTORSION,///<An improper torsion angle
};

///These are the recognized types of models
enum class Model_t{
    HARMONICOSCILLATOR,///<Bond-stretching treated harmonically
    FOURIERSERIES,///<A fourier series potential
    ELECTROSTATICS,///<A charge-charge interaction
    LENNARD_JONES,///<A 6-12 Lennard-Jones potential
};

///These are the recognized combination rules
enum class CombRule_t {
    PRODUCT,///<The product of the parameters \f$\prod_{i=1}^Nx_i\f$
    ARITHMETIC,///<A normal average \f$\frac{1}{N}\sum_{i=1}^Nx_i\f$
    GEOMETRIC,///<Geometric average \f$\left(\prod_{i=1}^Nx_i\right)^{1/N}\f$
};

///Flag for using the atom type or the atom class
enum class TypeTypes_t{
    TYPE,///<The term uses the atom type
    CLASS,///<The term uses the atom class
};

///Type of a FFTerm
using FFTerm_t=std::pair<Model_t,IntCoord_t>;

///Type of a quantity we are treating as a mathematical vector
using Vector=std::vector<double>;

///Shared pointer to a mathematical vector
using SharedVector=std::shared_ptr<Vector>;

///Shared pointer to a const mathematical vector
using cSharedVector=std::shared_ptr<const Vector>;

///Array of unsigned long integers
using IVector=std::vector<size_t>;

}//end namespace

//Instatiate some templates we know we're going to use
template class std::vector<double>;
template class std::vector<size_t>;
