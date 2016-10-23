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
/** \file Distance.hpp
 * 
 * \version 0.1
 * \date October 9, 2016 at 12:51 PM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FMANII_DISTANCE_HPP
#define FMANII_DISTANCE_HPP
#include "ForceManII/InternalCoordinates.hpp"

///Namespace for all code associated with ForceManII
namespace FManII {

///Implements derivatives for distance between points.  
///See [Distance Class](@ref distance) for more detail.
class Distance: public IntCoords {
public:
    Distance(std::shared_ptr<const VDouble> Carts):
        IntCoords(2,Carts){}
protected:
    VDouble compute_value(size_t deriv_i,Atoms_t coord_i)const;
};


} //End namespace FManII

#endif /* End header guard */

