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

namespace FManII {

inline Vector angle(const double* r1, const double* r2, const double* r3){
   const std::array<double,3> r12=diff(r1,r2),r32=diff(r3,r2);
   return {angle(r12,r32,cross(r12,r32))};
}


Vector Angle::deriv(size_t deriv_i,const Vector& sys,const IVector& coord_i)const{
    CHECK(deriv_i<1,"Higher order derivatives are not yet implemented!!!");
    const size_t atomi=coord_i[0],atomj=coord_i[1],atomk=coord_i[2];
    const double *q1=&(sys[atomi*3]),
                 *q2=&(sys[atomj*3]),
                 *q3=&(sys[atomk*3]);
    if(deriv_i==0) return angle(q1,q2,q3);
}

} //End namespace FManII

