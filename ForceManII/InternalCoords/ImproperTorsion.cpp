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
#include "ForceManII/InternalCoords/ImproperTorsion.hpp"
#include "ForceManII/Util.hpp"
#include "ForceManII/Common.hpp"
#include <algorithm>//For std::rotate

namespace FManII {


Vector ImproperTorsion::compute_value_(size_t deriv_i,const IVector& coord_i)const{
    std::array<size_t,3> atoms={coord_i[0],coord_i[2],coord_i[3]};
    size_t NDims=(size_t)std::pow(12.0,deriv_i);
    Vector phi(NDims,0.0);
    do{
        Vector temp=
        Torsion::compute_value_(deriv_i,{atoms[0],coord_i[1],atoms[1],atoms[2]});
        for(size_t i=0;i<NDims;++i)phi[i]+=(1.0/3.0)*temp[i];
        std::rotate(atoms.begin(),atoms.begin()+1,atoms.end());
    }while(*atoms.begin()!=coord_i[0]);
    return phi;
}

} //End namespace FManII
