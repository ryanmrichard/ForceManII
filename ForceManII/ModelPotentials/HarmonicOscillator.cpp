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

#include "ForceManII/ModelPotentials/HarmonicOscillator.hpp"
#include "ForceManII/Common.hpp"

namespace FManII{

Vector HarmonicOscillator::deriv(size_t order,
                                 const ParamInput_t& in_params,
                                 const CoordInput_t &in_coords)const
{
    const Vector &Qs=in_coords[0],&ks=in_params.at(Param_t::K),
                 &r0s=in_params.at(Param_t::r0);
    DEBUG_CHECK(ks.size()==r0s.size(),"len(Ks) != len(R0s)");
    DEBUG_CHECK(ks.size()==Qs.size(),"len(params) != len(coords)");
    const size_t N=Qs.size();
    const size_t NElems=(size_t)std::pow(N,order);
    Vector return_value(NElems,0.0);//Is correct already for Order>2
    if(order==0)
        for(size_t i=0;i<N;++i)
            return_value[0]+=0.5*ks[i]*(Qs[i]-r0s[i])*(Qs[i]-r0s[i]);
    else if(order==1)
        for(size_t i=0;i<N;++i)
            return_value[i]=ks[i]*(Qs[i]-r0s[i]);
    else if(order==2)
        for(size_t i=0;i<N;++i)
            return_value[i*N+i]=ks[i];
    return return_value;
}

}//End namespace
