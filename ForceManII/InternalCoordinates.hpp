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
#include "ForceManII/FManIIDefs.hpp"

namespace FManII {

class InternalCoordinates{
public:
    void add_coord(const IVector& atoms){
        atoms_.push_back(atoms);
        coords_->push_back(compute_value_(0,atoms)[0]);
    }

    const Vector& get_coords()const{return *coords_;}
    const std::vector<IVector>& get_types()const{return atoms_;}

    ///The name of this internal coordinate
    const std::string name;

protected:
    ///A vector such that element i is the value of the i-th coord
    SharedVector coords_;

    ///The Cartesian coordinates of the system
    cSharedVector carts_;

    ///A list such that element i is the NAtoms associated with the i-th coord
    std::vector<IVector> atoms_;

    InternalCoordinates(cSharedVector system,const std::string& namein):
        name(namein),coords_(std::make_shared<Vector>()),carts_(system){}

    ///Override in derived classes so that it returns the deriv of the coord
    virtual Vector compute_value_(size_t,const IVector&)const=0;


};

}//End namespace FManII

