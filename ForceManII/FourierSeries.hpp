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
/** \file FourierSeries.hpp
 * 
 * \version 0.1
 * \date November 19, 2016 at 2:09 PM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FMANII_FOURIERSERIES_HPP
#define FMANII_FOURIERSERIES_HPP

#include <vector>

///Namespace for all code associated with ForceManII
namespace FManII {

/** \brief Implements a Fourier Series model potential
 *
 * See [Fourier Series](@ref fourier) for details regarding conventions etc.
 *  
 */
struct FourierSeries {
    /** Computes the derivative of the energy of a FourierSeries
     * 
     * 
     *  \note This function expects all quantities to be in atomic units
     * 
     *  \param[in] Qs An N element vector such that the i-th element is the
     *                value of the i-th coordinate (torsion, impropertorsion,
     *                etc.) with the periodicity and phase shift already
     *                evalutate, i.e. give me n*phi-gamma for each angle
     *  \param[in] Vs An N element vector such that the i-th element is the
     *                value of the i-th coordinate's amplitude
     *  \param[in] ns A N element vector such that the i-th element is the value
     *                of the i-th coordinate's periodicity
     *  \param[in] Order What order derivative are we returning?
     *  \return The derivative in atomic units. 
     *
     */
    std::vector<double> deriv(const std::vector<double>& Qs,
                              const std::vector<double>& Vs,
                              const std::vector<double>& ns,
                              unsigned int Order=0)const;
};


} //End namespace FManII

#endif /* End header guard */

