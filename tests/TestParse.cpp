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

#include "ubiquitin.hpp"
#include <ForceManII/InternalCoordinates.hpp>
#include <iostream>
#include "TestMacros.hpp"


int main(int argc, char** argv){
    std::cout<<"Testing Input Parsing"<<std::endl;
 
    std::vector<double> FakeCarts({1.0});
    bool threw=false;
    try{
        FManII::CoordArray fakecoords=FManII::get_coords(FakeCarts,ubiquitin_FF_types,
                ubiquitin_FF_params,ubiquitin_conns);   
    }
    catch(...) {threw=true;}
    if(!threw)throw std::runtime_error("FManII internal check failed");

    
    FManII::CoordArray coords=FManII::get_coords(ubiquitin,ubiquitin_FF_types,
            ubiquitin_FF_params,ubiquitin_conns);
    
    const std::vector<double> &bonds=coords[FManII::Bond]->values();
    compare_vectors(bonds,ubiquitin_bonds,0.001,"generated bonds are not correct");
    
    const std::vector<double> &bondk=coords[FManII::Bond]->params(FManII::K);
    const std::vector<double> &bondr0=coords[FManII::Bond]->params(FManII::r0);
    
    compare_vectors(bondk,ubiquitin_K,0.001,"assigned Ks are not correct");
    compare_vectors(bondr0,ubiquitin_r0,0.001,"assigned r0s are not correct");
    
    
    return 0;
}

