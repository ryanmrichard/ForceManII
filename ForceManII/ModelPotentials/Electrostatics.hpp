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

#include "ForceManII/ModelPotential.hpp"

///Namespace for all code associated with ForceManII
namespace FManII {

/** \brief Implements the electrostatics model potential
 *
 * See [Electrostatics](@ref electrostatics) for details regarding conventions etc.
 *  
 */
struct Electrostatics: public ModelPotential {

    Electrostatics():ModelPotential({Param_t::q}){}

    /** Computes the derivative of the energy of according to Coulomb's law
     * 
     *  \note This function expects all quantities to be in atomic units
     * 
     *  \param[in] Order What order derivative are we returning?
     *  \param[in] in_params the qs
     *  \param[in] in_coords the value of the internal coordinates
     *  \return The derivative in atomic units. 
     *
     */
    Vector deriv(size_t order,
                 const ParamInput_t& in_params,
                 const CoordInput_t& in_coords)const;
};


} //End namespace FManII

