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
#include <ForceManII/ModelPotentials/HarmonicOscillator.hpp>
#include "TestMacros.hpp"
#include "ubiquitin.hpp"
#include <cmath>

using namespace FManII;
using namespace std;

int main(int argc, char** argv){

    test_header("Testing Harmonic Oscillator force-field term");
    HarmonicOscillator HO;
#ifndef NDEBUG
std::vector<double> a(2,0.0),d;
std::map<FManII::Param_t,std::vector<double>> ps={
    {FManII::Param_t::K,{0.0}},
    {FManII::Param_t::r0,{}}};
TEST_THROW(d=HO.deriv(0,ps,{a}),"len(k)!=len(r0)");
ps[FManII::Param_t::r0].push_back(0.0);
TEST_THROW(d=HO.deriv(0,ps,{a}),"len(k)!=len(r)");
#endif

    FManII::DerivType deriv=
            FManII::run_forcemanii(0,ubiquitin,ubiquitin_conns,
                                       FManII::amber99,ubiquitin_FF_types);
     for(auto qi : {IntCoord_t::BOND,IntCoord_t::ANGLE}){
        auto term_type=std::make_pair(FManII::Model_t::HARMONICOSCILLATOR,qi);
        //Energy check
        if(qi==IntCoord_t::BOND)
            test_value(deriv.at(term_type)[0],ubiquitinbond_e,1e-5,"Bond Energy");
        else
            test_value(deriv.at(term_type)[0],ubiquitinangle_e,1e-5,"Angle Energy");
    }
    test_footer();
    return 0;
} //End main

