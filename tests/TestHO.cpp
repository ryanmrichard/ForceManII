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
#include <ForceManII/FManIIDefs.hpp>
#include <ForceManII/ModelPotentials/HarmonicOscillator.hpp>
#include "TestMacros.hpp"
#include <cmath>

using namespace FManII;
using namespace std;

int main(int argc, char** argv){
    test_header("Testing Harmonic Oscillator force-field term");
    HarmonicOscillator HO;
    vector<double> a({3.0,2.0,1.0}),d;
    map<string,vector<double>> ps={{Param_t::K,{2.0}},{Param_t::r0,{}}};
#ifndef NDEBUG
TEST_THROW(d=HO.deriv(0,ps,{a}),"len(k)!=len(r0)");
ps[Param_t::r0].push_back(0.0);
TEST_THROW(d=HO.deriv(0,ps,{a}),"len(k)!=len(r)");
#endif
    ps[Param_t::r0]=vector<double>({1.5,1.0,0.5});
    ps[Param_t::K]=vector<double>({6.0,7.0,8.0});

    auto deriv=HO.deriv(0,ps,{a});
    test_value(deriv[0],11.25,1e-5,"Harmonic oscillator energy");
    deriv=HO.deriv(1,ps,{a});
    compare_vectors(deriv,{9.0,7.0,4.0},1e-5,"Gradient of harmonic oscillator");
    deriv=HO.deriv(2,ps,{a});
    compare_vectors(deriv,{6.0,0.0,0.0,0.0,7.0,0.0,0.0,0.0,8.0},
                    1e-5,"Hessian of harmonic oscillator");

    test_footer();
    return 0;
} //End main

