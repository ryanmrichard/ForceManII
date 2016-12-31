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
#include "ForceManII/InternalCoords/Torsion.hpp"
#include "ForceManII/Util.hpp"
#include "ForceManII/Common.hpp"

using DArray=std::array<double,3>;
namespace FManII {

inline Vector torsion(const double* q1, const double* q2,
        const double* q3, const double* q4){
    const DArray r21=diff(q1,q2),r32=diff(q2,q3),r34=diff(q4,q3);
    const DArray n1=cross(r21,r32),n2=cross(r34,r32);
    return {angle(n1,n2,r32)};
}

Vector Torsion::compute_value_(size_t deriv_i,const IVector& coord_i)const{
    CHECK(deriv_i<1,"Higher order derivatives are not yet implemented!!!");
    const size_t atomi=coord_i[0],atomj=coord_i[1],
                 atomk=coord_i[2],atoml=coord_i[3];
    const double *q1=&((*carts_)[atomi*3]), *q2=&((*carts_)[atomj*3]),
                 *q3=&((*carts_)[atomk*3]), *q4=&((*carts_)[atoml*3]);
    if(deriv_i==0) return torsion(q1,q2,q3,q4);
}


} //End namespace FManII

