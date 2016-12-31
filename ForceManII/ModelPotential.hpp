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
#include "ForceManII/InternalCoordinates.hpp"
#include <vector>
#include <array>

namespace FManII {


struct ModelPotential{
    using ParamInput_t=std::map<Param_t,Vector>;
    using CoordInput_t=std::vector<Vector>;
    ModelPotential(const std::vector<Param_t>& ps):
        params(ps){}

    ///The types of the parameters
    const std::vector<Param_t> params;

    bool operator==(const ModelPotential& other)const{
        return params==other.params;
    }

    bool operator!=(const ModelPotential& other)const{
        return !this->operator==(other);
    }

    ///The function to override so it returns your derivative
    virtual Vector deriv(size_t order,
                         const ParamInput_t& in_params,
                         const CoordInput_t& in_coords)const=0;

};

}//End namespace
