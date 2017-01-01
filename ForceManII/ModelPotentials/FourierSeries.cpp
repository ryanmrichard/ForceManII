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

#include "ForceManII/ModelPotentials/FourierSeries.hpp"
#include "ForceManII/Common.hpp"
#include <cmath>

namespace FManII{

namespace detail{
inline double odd_deriv_(size_t d,double Q,double Q0, double V, double n){
    const double de=(n==d && n==0.0? 0.0: V*std::pow(n,d)*std::sin(n*Q-Q0));
    return (d+1/2%2==0?de:-1*de);
}

inline double even_deriv_(size_t d,double Q, double Q0, double V, double n){
    const double de=(n==d && n==0.0? 0.0 : V*std::pow(n,d)*std::cos(n*Q-Q0));
    return (d/2%2==0?de:-1*de);
}

}//End detail

Vector FourierSeries::deriv(size_t order,
                            const ParamInput_t &in_params,
                            const CoordInput_t &in_coords)const
{
    DEBUG_CHECK(order<3,"Despite knowing the form for derivs higher than 2,"
                        " I haven't coded them up");
    const Vector &Qs=in_coords[0],&Vs=in_params.at(Param_t::amp),
                 &phis=in_params.at(Param_t::phi),
                 &ns=in_params.at(Param_t::n);
    const size_t N=Qs.size(),Nps=Vs.size();
    DEBUG_CHECK(Nps==phis.size(),"Number of amplitudes != number of phis");
    DEBUG_CHECK(Nps==ns.size(),"Number of amplitudes != number of ns");
    DEBUG_CHECK(Nps%N==0,"Number of parameters not an integer multiple of the"
                         "number of coordinates");
    const size_t dim=Nps/N;
    const size_t NElems=(size_t)std::pow(N,order);
    Vector return_value(NElems,0.0);
    for(size_t i=0;i<N;++i){
        for(size_t j=0;j<dim;++j){
            const size_t idx=i*dim+j;
            if(order==0)
                return_value[0]+=Vs[idx]+
                        detail::even_deriv_(0,Qs[i],phis[idx],Vs[idx],ns[idx]);
            else if(order==1)
                return_value[i]=
                        detail::odd_deriv_(1,Qs[i],phis[idx],Vs[idx],ns[idx]);
            else if(order==2){
                return_value[i*N+i]=
                        detail::even_deriv_(2,Qs[i],phis[idx],Vs[idx],ns[idx]);
                for(size_t k=0;k<i;++k)
                    for(size_t l=0;l<dim;++l){
                        const size_t idx2=k*dim+l;
                        return_value[i*N+k]=return_value[k*N+i]=
                                detail::odd_deriv_(1,Qs[i],phis[idx],Vs[idx],ns[idx])+
                                detail::odd_deriv_(1,Qs[k],phis[idx2],Vs[idx2],ns[idx2]);
                    }
            }
        }

    }
    return return_value;
}
}//End namespace
