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

#include "ForceManII/ModelPotentials/Electrostatics.hpp"
#include "ForceManII/Common.hpp"
#include <cmath>
namespace FManII{

Vector Electrostatics::deriv(size_t order,
                             const ParamInput_t& in_params,
                             const CoordInput_t &in_coords)const
{
    const Vector &Qs=in_coords[0],&qs=in_params.at(Param_t::q);
    const size_t n=Qs.size();
    DEBUG_CHECK(n==qs.size(),"Qs must be same length as qs");
    DEBUG_CHECK(order<=2,"Derivatives larger than order 2 are not coded");

    Vector d(static_cast<size_t>(std::pow(n,order)),0.0);
    if(order==0)
        for(size_t i=0;i<n;++i)d[0]+=qs[i]/Qs[i];
    else if(order==1)
        for(size_t i=0;i<n;++i)d[i]=-qs[i]/(Qs[i]*Qs[i]);
    else if(order==2)
        for(size_t i=0;i<n;++i)d[i*n+i]=2.0*qs[i]/(Qs[i]*Qs[i]*Qs[i]);
    return d;
}

}//End namespace
