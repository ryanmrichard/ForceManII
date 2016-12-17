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
#include "ForceManII/Distance.hpp"
#include "ForceManII/Common.hpp"
#include "ForceManII/Util.hpp"

namespace FManII {

inline std::vector<double> dist(const double* q1, const double* q2,double off){
    return {mag(diff(q1,q2))-off};
}

Distance::VDouble Distance::compute_value(size_t deriv_i,Atoms_t coord_i)const{
    CHECK(deriv_i<1,"Higher order derivatives are not yet implemented!!!");
    const size_t atomi=coord_i[0],atomj=coord_i[1];
    const double* q1=&((*carts_)[atomi*3]), *q2=&((*carts_)[atomj*3]);
    const double off=params_.count(r0) ? params_.at(r0).back() : 0.0;
    if(deriv_i==0) return dist(q1,q2,off);
}

} //End namespace FManII

