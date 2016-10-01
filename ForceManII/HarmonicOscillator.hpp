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
/** \file HarmonicOscillator.hpp
 * 
 * \version 0.1
 * \date October 1, 2016 at 12:33 PM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */
#ifndef FMANII_HARMONICOSCILLATOR_HPP
#define FMANII_HARMONICOSCILLATOR_HPP

#include <vector>

///Namespace for all code associated with ForceManII
namespace FManII {

/** Implements the harmonic oscillator model potential
 *
 *  The harmonic oscillator (HO) is one of the most common approximate
 *  potentials in all of physics.  The energy of a \f$N\f$-dimensional 
 *  harmonic oscillator is given by:
 * 
 *  \f[
 *      E(\mathbf{q};\mathbf{k})=\frac{1}{2}\mathbf{k}\dot 
 *            \left(\mathbf{q}\circ\mathbf{q}\right),
 *  \f]
 *  where \f$\mathbf{q}\f$ is a vector of \f$N\f$ coordinates and 
 *  \f$\mathbf{k}\f$ is a vector of
 *  \f$N\f$ force constants such that the \f$q_i\f$ has force constant 
 *  \f$k_i\f$.  \f$\circ\f$ denotes the Hadamard product (element-wise
 *  multiplication), i.e. \f$\left[\mathbf{q}\circ\mathbf{q}\right]_{i}=q_i^2\f$.
 * 
 *  The derivative w.r.t. to \f$\mathbf{q}\f$ is:
 *  \f[
 *      \frac{\partial E(\mathbf{q};\mathbf{k})}{\partial \mathbf q}=
 *       \mathbf{k}\dot \mathbf{q},
 *  \f]
 *  the Hessian is just \f$\mathbf{k}\f$ down the diagonal.  All higher
 *  derivatives are 0.  Doesn't get much easier than that.
 * 
 */

struct HarmonicOscillator {
    /** Computes the geometric derivative of the energy of a Harmonic Oscillator
     * 
     *  \note This function expects all quantities to be in atomic units
     * 
     *  \param[in] Qs An N element vector such that the i-th element is the
     *                value of the i-th coordinate (bond-length, angle, etc.)
     *  \param[in] Ks An N element vector such that the i-th element is the
     *                value of the i-th coordinate's force constant
     *  \param[in] Order What order derivative are we returning?
     *  \return The derivative in atomic units. 
     *
     */
    std::vector<double> deriv(const std::vector<double>& Qs,
                              const std::vector<double>& Ks,
                              unsigned int Order=0)const;
};

} //End namespace FManII

#endif /* End header guard */

