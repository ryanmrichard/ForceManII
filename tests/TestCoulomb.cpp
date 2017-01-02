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
#include <ForceManII/FManII.hpp>
#include <ForceManII/ModelPotentials/Electrostatics.hpp>
#include "TestMacros.hpp"
#include "testdata/ubiquitin.hpp"
#include <cmath>

int main(int argc, char** argv){
    test_header("Testing charge-charge force-field term");
    FManII::Electrostatics ct;
    
#ifndef NDEBUG
std::vector<double> a(3,0.0),b(2,0.0),c;
TEST_THROW(c=ct.deriv(0,{{FManII::Param_t::q,a}},{b}),
           "charges.size()==dist.size()");
#endif

    FManII::DerivType deriv=FManII::run_forcemanii(0,
                  ubiquitin,ubiquitin_conns,FManII::amber99,ubiquitin_FF_types);
    for(auto param:{FManII::IntCoord_t::PAIR,FManII::IntCoord_t::PAIR14}){
         const FManII::FFTerm_t term_type=
                std::make_pair(FManII::Model_t::ELECTROSTATICS,param);
        //There appears to be a slight loss in precision for this term...
        test_value(deriv.at(term_type)[0],ubiquitin_egys.at(term_type),4e-4,"Charge-Charge Energy");
     }

        
    test_footer();
    return 0;
} //End main
