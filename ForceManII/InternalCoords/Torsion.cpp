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
    const DArray r21=diff(q2,q1),r23=diff(q2,q3),r34=diff(q3,q4);
    const DArray n1=cross(r21,r23),n2=cross(r34,r23);
    return {angle(n1,n2)};
}

//Gradient of the dot product
inline Vector ddot(const double* q1,const double* q2,
                   const double* q3,const double* q4){
    const DArray r21=diff(q2,q1),r23=diff(q2,q3),r34=diff(q3,q4);
    const DArray r31=diff(q3,q1),r42=diff(q4,q2);
    const DArray n1=cross(r21,r23),n2=cross(r34,r23);
    const DArray ddr1=cross(n2,r23),ddr2=sum(cross(n1,r34),cross(n2,r31));
    const DArray ddr3=diff(cross(n1,r42),cross(n2,r21)),ddr4=cross(n1,r23);
    return {ddr1[0],ddr1[1],ddr1[2],ddr2[0],ddr2[1],ddr2[2],
            ddr3[0],ddr3[1],ddr3[2],ddr4[0],ddr4[1],ddr4[2]};
}

//Gradient of A dot the cross product
inline Vector dcross(const double* q1,const double* q2,
                     const double* q3,const double* q4)
{
    const DArray r21=diff(q2,q1),r23=diff(q2,q3),r34=diff(q3,q4);
    const DArray r31=diff(q3,q1),r42=diff(q4,q2);
    const DArray n1=cross(r21,r23),n2=cross(r34,r23);
    const DArray A=cross(n1,n2);
    const double r23dA=dot(r23,A),r23dn2=dot(r23,n2),r21dA=dot(r21,A);
    const double r31dA=dot(r31,A),r34dA=dot(r34,A),r42dA=dot(r42,A);
    const double r23dn1=dot(r23,n1);
    const double x=dot(r31,n2)-dot(r34,n1),y=dot(r21,n2)+dot(r42,n1);

    const DArray dcr1{A[0]*r23dn2-n2[0]*r23dA,
                      A[1]*r23dn2-n2[1]*r23dA,
                      A[2]*r23dn2-n2[2]*r23dA};
    const DArray dcr2{A[0]*x-n2[0]*r31dA+n1[0]*r34dA,
                      A[1]*x-n2[1]*r31dA+n1[1]*r34dA,
                      A[2]*x-n2[2]*r31dA+n1[2]*r34dA
    };
    const DArray dcr3{-A[0]*y+n2[0]*r21dA+n1[0]*r42dA,
                      -A[1]*y+n2[1]*r21dA+n1[1]*r42dA,
                      -A[2]*y+n2[2]*r21dA+n1[2]*r42dA
    };
    const DArray dcr4{-A[0]*r23dn1+n1[0]*r23dA,
                      -A[1]*r23dn1+n1[1]*r23dA,
                      -A[2]*r23dn1+n1[2]*r23dA};
    return  {dcr1[0],dcr1[1],dcr1[2],dcr2[0],dcr2[1],dcr2[2],
             dcr3[0],dcr3[1],dcr3[2],dcr4[0],dcr4[1],dcr4[2]};
}

inline Vector dtorsion(const double* q1,const double* q2,
                       const double* q3,const double* q4)
{
    const DArray r21=diff(q2,q1),r23=diff(q2,q3),r34=diff(q3,q4);
    const DArray n1=cross(r21,r23),n2=cross(r34,r23);
    const DArray A=cross(n1,n2);
    const double n1dn2=dot(n1,n2),magA=mag(A);
    const double pf=1.0/(magA*magA+n1dn2*n1dn2),tanphi=magA/n1dn2;
    const Vector dot_part=ddot(q1,q2,q3,q4);
    const Vector c_part=dcross(q1,q2,q3,q4);
    Vector rv(12);
    for(size_t i=0;i<12;++i)rv[i]=pf*(c_part[i]/tanphi-magA*dot_part[i]);
    return rv;
}


Vector Torsion::deriv(size_t deriv_i,const Vector& sys,const IVector& coord_i)const{
    CHECK(deriv_i<2,"Higher order derivatives are not yet implemented!!!");
    const size_t atomi=coord_i[0],atomj=coord_i[1],
                 atomk=coord_i[2],atoml=coord_i[3];
    const double *q1=&(sys[atomi*3]), *q2=&(sys[atomj*3]),
                 *q3=&(sys[atomk*3]), *q4=&(sys[atoml*3]);
    if(deriv_i==0) return torsion(q1,q2,q3,q4);
    if(deriv_i==1) return dtorsion(q1,q2,q3,q4);
}


} //End namespace FManII

