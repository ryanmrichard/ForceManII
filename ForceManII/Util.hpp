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
/** \file Util.hpp
 * 
 * Contains helper utility functions of a mathematical nature
 * 
 * \version 0.1
 * \date November 6, 2016 at 1:38 PM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FMANII_UTIL_HPP
#define FMANII_UTIL_HPP
#include <array>
#include <numeric>
#include <cmath>

///Namespace for all code associated with ForceManII
namespace FManII {

///Returns the cross product of two vectors
template<typename T>
std::array<double,3> cross(const T& v1, const T& v2){
    return {v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0]};
}
///Returns the difference between two vectors, v1-v2
template<typename T>
std::array<double,3> diff(const T& v1,const T& v2){
    return {v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]};
}


///Returns the dot product of two vectors
inline double dot(const std::array<double,3>& v1,
                  const std::array<double,3>& v2){
    return std::inner_product(v1.begin(),v1.end(),v2.begin(),0.0);
}

///Returns the magnitude of a vector
inline double mag(const std::array<double,3>& v1){
    return std::sqrt(dot(v1,v1));
}

///Returns the angle between two vectors given the vectors \p v1, \p v2 and a
///normal to them \p n (vectors need not be units already)
inline double angle(const std::array<double,3>& v1,
                    const std::array<double,3>& v2,
                    const std::array<double,3>& n){
    const double mag1=mag(v1),mag2=mag(v2),magn=mag(n);
    const std::array<double,3> v1crossv2=cross(v1,v2);
    const double v1dotv2=dot(v1,v2),v1cv2dotn=dot(v1crossv2,n);
    const double cosphi=v1dotv2/(mag1*mag2),sinphi=v1cv2dotn/(mag1*mag2*magn);
    return std::atan2(sinphi,cosphi);   
}


} //End namespace FManII

#endif /* End header guard */

