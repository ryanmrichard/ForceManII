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

#include "ForceManII/InternalCoordinates.hpp"

///Namespace for all code associated with ForceManII
namespace FManII {

///Implements derivatives for distance between points.  
///See [Distance Class](@ref distance) for more detail.
class Distance: public InternalCoordinates {
public:
    Distance(const std::string& namein):
        InternalCoordinates(namein){}
    Vector deriv(size_t deriv_i,const Vector& sys,const IVector& coord_i)const;
};

class Bond:public Distance{
public:
    Bond():Distance(IntCoord_t::BOND){}
};

class Pair:public Distance{
public:
    Pair():Distance(IntCoord_t::PAIR){}
};

///We assume that this is being used in a Urey-Bradley term *i.e* that the
///vertex of the angle is being given to us as well
class Pair13: public Distance{
public:
    Pair13():Distance(IntCoord_t::PAIR13){}
    Vector deriv(size_t deriv_i,const Vector& sys,const IVector& coord_i)const;
};

class Pair14: public Distance{
public:
    Pair14():Distance(IntCoord_t::PAIR14){}
};

} //End namespace FManII


