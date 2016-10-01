/*
 * Copyright (C) 2016 Ryan M. Richard.
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
/** \file FFTerm.hpp
 * 
 * \version 0.1
 * \date October 1, 2016 at 11:32 AM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FFTERM_HPP
#define FFTERM_HPP

///Namespace for all code associated with ForceManII
namespace FManII {


/** The class responsible for computing derivatives of a force field term
 * 
 *  All energetic terms in a forcefield depend on some generic coordinate
 *  \f$q\f$, which itself depends on the standard Cartesian coordinates of the
 *  system.  Consequentially, the first derivative of the energy with respect to
 *  one of these coordinates, say \f$x\f$ is:
 *  \f[
      \frac{\partial E(q)}{\partial x}=\frac{\partial f(q)}{\partial q}
                                         \frac{\partial q}{\partial x},
    \f]
 *  where \f$f(q)\f$ is the functional form of the term (harmonic oscillator, 
 *  Fourier series, etc).  The second derivative with respect to \f$x\f$ and
 *  \f$y\f$ is then:
 *  \f[
      \frac{\partial^2 E(q)}{\partial x\partial y}=
        \frac{\partial^2 f(q)}{\partial q^2}\frac{\partial q}{\partial x}
        \frac{\partial q}{\partial y}+
        \frac{\partial f(q)}{\partial q}
        \frac{\partial q^2}{\partial x\partial y}
    \f]
 *  
 *  This class is responsible for handeling these derivative manipulations in a
 *  generic way.  
 * 
 */
template <typename Model_t, typename Coord_t>
class FFTerm {
};

} //End namespace FManII

#endif /* End header guard */

