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
#include "ForceManII/LennardJones.hpp"
#include "ForceManII/Common.hpp"
namespace FManII {

std::vector<double> LennardJones::deriv(const std::vector<double>& Qs,
                                        const std::vector<double>& sigmas,
                                        const std::vector<double>& epsilons,
                                        unsigned int Order)const{
    DEBUG_CHECK(Qs.size()==sigmas.size(),"Qs must be same length as sigmas");
    DEBUG_CHECK(Qs.size()==epsilons.size(),"Qs must be same length as epsilons");

    std::vector<double> d(1,0.0);
    if(Order==0)
        for(size_t i=0;i<Qs.size();++i){
            const double term=sigmas[i]/Qs[i];
            const double term2=term*term;
            const double term6=term2*term2*term2;
            d[0]+=epsilons[i]*(term6*term6-2.0*term6);
        }
    return d;
}

} //End namespace FManII
