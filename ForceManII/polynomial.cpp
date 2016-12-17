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

namespace FManII {

Return_t Polynomial::deriv(const Return_t& Qs, 
                                   const Return_t& Ks, 
                                   unsigned int Order) const {
    DEBUG_CHECK(Qs.size()==Ks.size(),"Qs must be same length as Ks");
    const size_t N=Qs.size();
    const size_t NElems=(size_t)std::pow(N,Order);
    Return_t return_value(NElems,0.0);//Is correct already for Order>2
    if(Order==0)
        for(size_t i=0;i<N;++i)return_value[0]+=0.5*Ks[i]*Qs[i]*Qs[i];
    else if(Order==1)
        for(size_t i=0;i<N;++i)return_value[i]=Ks[i]*Qs[i];
    else if(Order==2)
        for(size_t i=0;i<N;++i)return_value[i*N+i]=Ks[i];
    return return_value;
}


} //End namespace FManII

