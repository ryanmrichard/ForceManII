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

#include "ForceManII/InternalCoords/Angle.hpp"
#include "ForceManII/Common.hpp"
#include "ForceManII/Util.hpp"
#include <iostream>

using namespace std;

namespace FManII {

inline Vector angle(const double* r1, const double* r2, const double* r3)
{
   const array<double,3> r12=diff(r1,r2),r32=diff(r3,r2);
   return {angle(r12,r32)};
}

inline Vector dangle(const double* r1,const double* r2, const double* r3)
{
    const array<double,3> r12=diff(r1,r2),r31=diff(r3,r1),r32=diff(r3,r2);
    const array<double,3> n=cross(r12,r32);
    const double r12_d_r32=dot(r12,r32);
    const double magn=mag(n);
    const double tantheta=r12_d_r32/magn;//Actually 1 over tan(theta)
    const array<double,3> A=cross(r32,n),B=cross(n,r31),C=cross(n,r12);
    const double prefactor=1/(std::pow(magn,2)+std::pow(r12_d_r32,2));
    return {prefactor*(tantheta*A[0]-r32[0]*magn),
            prefactor*(tantheta*A[1]-r32[1]*magn),
            prefactor*(tantheta*A[2]-r32[2]*magn),
            prefactor*(tantheta*B[0]+(r12[0]+r32[0])*magn),
            prefactor*(tantheta*B[1]+(r12[1]+r32[1])*magn),
            prefactor*(tantheta*B[2]+(r12[2]+r32[2])*magn),
            prefactor*(tantheta*C[0]-r12[0]*magn),
            prefactor*(tantheta*C[1]-r12[1]*magn),
            prefactor*(tantheta*C[2]-r12[2]*magn)
    };
}

Vector Angle::deriv(size_t deriv_i,const Vector& sys,const IVector& coord_i)const{
    CHECK(deriv_i<2,"Higher order derivatives are not yet implemented!!!");
    const size_t atomi=coord_i[0],atomj=coord_i[1],atomk=coord_i[2];
    const double *q1=&(sys[atomi*3]),
                 *q2=&(sys[atomj*3]),
                 *q3=&(sys[atomk*3]);
    if(deriv_i==0) return angle(q1,q2,q3);
    if(deriv_i==1) return dangle(q1,q2,q3);
}

} //End namespace FManII

