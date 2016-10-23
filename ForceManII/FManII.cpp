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
#include <iostream>

namespace FManII {
using DVector=std::vector<double>;
using IVector=std::vector<size_t>;
template<size_t n> using IArray=std::array<size_t,n>;
using shared_DVector=std::shared_ptr<const DVector>;
using std::max;
using std::min;

CoordArray get_coords(const DVector& Carts,const AtomTypes& Types,
                     const ParamTypes& Params,const ConnData& Conns){
    const size_t NAtoms=Carts.size()/3;
    DEBUG_CHECK(NAtoms==Conns.size(),"Number of atoms differs among inputs");
    shared_DVector Sys=std::make_shared<DVector>(Carts);
    
    CoordArray FoundCoords;
    FoundCoords.emplace(Bond,std::move(make_unique<Distance>(Sys)));
    
    for(size_t AtomI=0;AtomI<NAtoms;++AtomI){
        for(size_t AtomJ : Conns[AtomI]){
            DEBUG_CHECK(AtomI!=AtomJ,"AtomI is bonded to itself. What??");
            if(AtomJ<AtomI)continue;
            const IVector pair({AtomI,AtomJ});
            const IArray<2> ts={Types[AtomI][0],Types[AtomJ][0]};
            const IArray<2> sorted_types={min(ts[0],ts[1]),max(ts[0],ts[1])};
            const IVector ff_types({sorted_types[0],sorted_types[1]});
            const double K_in=Params.at(Bond).at(K).at(ff_types);
            const double r0_in=Params.at(Bond).at(r0).at(ff_types);
            FoundCoords[Bond]->add_coord(pair,Param_t::K,K_in,Param_t::r0,r0_in);
        }
        
    }
    return FoundCoords;
}


} //End namespace FManII

