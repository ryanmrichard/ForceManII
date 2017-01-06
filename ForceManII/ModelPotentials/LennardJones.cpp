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

#include "ForceManII/ModelPotentials/LennardJones.hpp"
#include "ForceManII/Common.hpp"
#include <cmath>

namespace FManII{

Vector LennardJones::deriv(size_t order,
              const ParamInput_t &in_params,
              const CoordInput_t &in_coords)const
{
    const Vector &Qs=in_coords[0],
                 &ss=in_params.at(FManII::Param_t::sigma),
                 &es=in_params.at(FManII::Param_t::epsilon);
    const size_t n=Qs.size();
    DEBUG_CHECK(ss.size()==n,"len(Qs) != len(sigmas)");
    DEBUG_CHECK(es.size()==n,"len(Qs) != len(epsilons)");
    std::vector<double> d(static_cast<size_t>(std::pow(n,order)),0.0);
    if(order==0)
        for(size_t i=0;i<n;++i){
            const double term=ss[i]/Qs[i];
            const double term2=term*term;
            const double term6=term2*term2*term2;
            d[0]+=es[i]*term6*(term6-2.0);
        }
    else if(order==1)
        for(size_t i=0;i<n;++i){
            const double term=ss[i]/Qs[i];
            const double term2=term*term;
            const double term6=term2*term2*term2;
            d[i]+=12.0*es[i]*term6/Qs[i]*(1.0-term6);
        }
    else if(order==2)
        for(size_t i=0;i<n;++i){
            const double term=ss[i]/Qs[i];
            const double term2=term*term;
            const double term6=term2*term2*term2;
            d[i*n+i]+=12.0*es[i]*term6/(Qs[i]*Qs[i])*(13.0*term6-7.0);
        }
    return d;
}

}//End namespace
