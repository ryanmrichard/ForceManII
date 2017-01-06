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
#include <ForceManII/ModelPotentials/FourierSeries.hpp>
#include "TestMacros.hpp"
#include "testdata/ubiquitin.hpp"
#include <cmath>
using namespace std;
using namespace FManII;

int main(int argc, char** argv){
    test_header("Testing FourierSeries force-field term");
    FourierSeries FS;
    Vector theta({2.0,3.10,0.6}),d;
    map<string,Vector> ps={
        {Param_t::amp,{0.0}},
        {Param_t::phi,{0.0}},
        {Param_t::n,{0.0}}};
#ifndef NDEBUG
TEST_THROW(d=FS.deriv(0,ps,{theta}),"thetas.size()!=amps.size()");
ps[Param_t::amp]=Vector({3.2,2.2,1.2});
TEST_THROW(d=FS.deriv(0,ps,{theta}),"theta.size()!=phis.size()");
ps[Param_t::phi]=Vector({6.0,7.0,8.0});
TEST_THROW(d=FS.deriv(0,ps,{theta}),"thetas.size()!=ns.size()");
#endif
    ps[Param_t::amp]=Vector({3.2,2.2,1.2});
    ps[Param_t::phi]=Vector({M_PI/2.0,M_PI,M_PI/4});
    ps[Param_t::n]=Vector({2.0,3.0,4.0});
    const Vector egy({6.308577929519}),
        grad({-4.18331917,0.8213992,-4.79539532}),
        hess({9.687071941,0.0,0.0,0.0,-19.64606144,0.0,0.0,0.0,0.84079682});
    
    compare_vectors(FS.deriv(0,ps,{theta}),egy,1e-5,"Fourier series energy");
    compare_vectors(FS.deriv(1,ps,{theta}),grad,1e-5,"Fourier series gradient");
    compare_vectors(FS.deriv(2,ps,{theta}),hess,1e-5,"Fourier series Hessian");

    test_footer();
    return 0;
} //End main


