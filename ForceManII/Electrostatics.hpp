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
/** \file Electrostatics.hpp
 * 
 * \version 0.1
 * \date November 27, 2016 at 2:10 PM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FMANII_ELECTROSTATICS_HPP
#define FMANII_ELECTROSTATICS_HPP

#include<vector>

///Namespace for all code associated with ForceManII
namespace FManII {

/** \brief Implements the electrostatics model potential
 *
 * See [Electrostatics](@ref electrostatics) for details regarding conventions etc.
 *  
 */
struct Electrostatics {
    /** Computes the derivative of the energy of according to Coulomb's law
     * 
     *  \note This function expects all quantities to be in atomic units
     * 
     *  \param[in] Qs An N element vector such that the i-th element is the
     *                value of the i-th coordinate (bond-length, angle, etc.)
     *  \param[in] qs An N element vector such that the i-th element is the
     *                (possibly scaled) product of charges for the i-th coordinate
     *  \param[in] Order What order derivative are we returning?
     *  \return The derivative in atomic units. 
     *
     */
    std::vector<double> deriv(const std::vector<double>& Qs,
                              const std::vector<double>& qs,
                              unsigned int Order=0)const;
};

} //End namespace FManII

#endif /* End header guard */

