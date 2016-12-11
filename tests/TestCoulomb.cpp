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
#include <ForceManII/FFTerm.hpp>
#include <ForceManII/Electrostatics.hpp>
#include "TestMacros.hpp"
#include "ubiquitin.hpp"
#include <cmath>

int main(int argc, char** argv){
    test_header("Testing charge-charge force-field term");
    FManII::Electrostatics ct;
    
#ifndef NDEBUG
std::vector<double> a(3,0.0),b(2,0.0),c;
TEST_THROW(c=ct.deriv(a,b),"charges.size()==dist.size()");
#endif

     FManII::CoordArray coords=FManII::get_coords(ubiquitin,ubiquitin_FF_types,
            ubiquitin_FF_params,ubiquitin_conns,1/1.2);
     const std::vector<double>& dist=
        coords[FManII::ELECTROSTATICS]->values();
     const std::vector<double>& charges=
        coords[FManII::ELECTROSTATICS]->params(FManII::q);
    
    //Energy check
    std::vector<double> Energy=ct.deriv(dist,charges);
    //There appears to be a slight loss in precision for this term...
    test_value(Energy[0],ubiquitinelectrostatics_e,2e-4,"Charge-Charge Energy");
        
    test_footer();
    return 0;
} //End main
