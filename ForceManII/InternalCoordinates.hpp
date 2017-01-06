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

struct InternalCoordinates{
    InternalCoordinates(const std::string& name_):
        name(name_){}

    virtual Vector deriv(size_t order,const Vector& sys,const IVector& atoms)const=0;

    ///The name of this internal coordinate
    const std::string name;
};
}//End namespace FManII

