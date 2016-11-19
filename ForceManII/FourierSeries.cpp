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
#include "ForceManII/FourierSeries.hpp"
#include "ForceManII/Common.hpp"//For debug check
#include <cmath>
#include <iostream>

using Return_t=std::vector<double>;
namespace FManII {

inline double odd_deriv_(size_t d,double Q,double V, double n){
    const double de=V*std::pow(n,d)*std::sin(Q);
    return (d+1/2%2==0?de:-1*de);
}

inline double even_deriv_(size_t d,double Q, double V, double n){
    const double de=V*std::pow(n,d)*std::cos(Q);
    return (d/2%2==0?de:-1*de);
}

Return_t FourierSeries::deriv(const Return_t& Qs, 
                              const Return_t& Vs, 
                              const Return_t& ns,
                              unsigned int Order) const {
    DEBUG_CHECK(Order<3,"Despite knowing the form for derivs higher than 2,"
                        " I haven't coded them up");
    DEBUG_CHECK(Qs.size()==Vs.size(),"Qs must be same length as Vs");
    DEBUG_CHECK(Qs.size()==ns.size(),"Qs must be same length as ns");
    const size_t N=Qs.size();
    const size_t NElems=(size_t)std::pow(N,Order);
    Return_t return_value(NElems,0.0);//Is correct already for Order>2
    if(Order==0)
        for(size_t i=0;i<N;++i)
            return_value[0]+=Vs[i]+even_deriv_(0,Qs[i],Vs[i],ns[i]);
    else if(Order==1)
        for(size_t i=0;i<N;++i)
            return_value[i]=odd_deriv_(1,Qs[i],Vs[i],ns[i]);
    else if(Order==2)
        for(size_t i=0;i<N;++i){
            return_value[i*N+i]=even_deriv_(2,Qs[i],Vs[i],ns[i]);
            for(size_t j=0;j<i;++j)
                return_value[i*N+j]=return_value[j*N+i]=
                    odd_deriv_(1,Qs[i],Vs[i],ns[i])+
                    odd_deriv_(1,Qs[j],Vs[j],ns[j]);
        }
    return return_value;
}

} //End namespace FManII

