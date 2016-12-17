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
/** \file LennardJones.hpp
 *
 * \version 0.1
 * \date December 10, 2016 at 11:29 AM (EST)
 *
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 *
 * Additional contributions by:
 *
 */

#ifndef FMANII_LENNARDJONES_HPP
#define FMANII_LENNARDJONES_HPP

#include<vector>

///Namespace for all code associated with ForceManII
namespace FManII {

/** \brief Implements the Lennard-Jones 6-12 model potential
 *
 * See [Lennard-Jones](@ref lj) for details regarding conventions etc.
 *
 */
struct LennardJones {
    /** Computes the derivative of the energy of according to the 6-12
     *  Lennard-Jones potential
     *
     *  \note This function expects all quantities to be in atomic units
     *
     *  \param[in] Qs An N element vector such that the i-th element is the
     *                value of the i-th coordinate (bond-length, angle, etc.)
     *  \param[in] sigmas An N element vector such that the i-th element is the
     *                averaged sigmas for the i-th coordinate
     *  \param[in] epsilons An N element vector such that the i-th element is the
     *                averaged epsilons for the i-th coordinate
     *  \param[in] Order What order derivative are we returning?
     *  \return The derivative in atomic units.
     *
     */
    std::vector<double> deriv(const std::vector<double>& Qs,
                              const std::vector<double>& sigmas,
                              const std::vector<double>& epsilons,
                              unsigned int Order=0)const;
};

} //End namespace FManII

#endif /* End header guard */
