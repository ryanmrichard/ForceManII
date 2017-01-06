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
#include <ForceManII/ModelPotentials/LennardJones.hpp>
#include "TestMacros.hpp"
#include "testdata/ubiquitin.hpp"


using namespace std;
using namespace FManII;

int main(int argc, char** argv){
    test_header("Testing 6-12 force-field term");
    LennardJones lj;
    map<string,Vector> ps={
        {Param_t::sigma,{0.0}},
        {Param_t::epsilon,{0.0}}};
    Vector qs({3.2,2.2,1.2}),c;
#ifndef NDEBUG
TEST_THROW(c=lj.deriv(0,ps,{qs}),"sigmas.size()==dist.size()");
ps[FManII::Param_t::sigma]=Vector({M_PI/2.0,M_PI,M_PI/4.0});
TEST_THROW(c=lj.deriv(0,ps,{qs}),"epsilons.size()==dist.size()");
#endif
    ps[FManII::Param_t::sigma]=Vector({M_PI/2.0,M_PI,M_PI/4.0});
    ps[FManII::Param_t::epsilon]=Vector({2.0,3.0,4.0});
    const Vector egy({164.162849653}),
                 grad({1.03457493e-01,-1.03778526e+03,2.89706016e+00}),
        hess({-2.23560931e-01,0.0,0.0,0.0,6.51078520e+03,0.0,0.0,0.0,-1.56637591e+01});
    compare_vectors(lj.deriv(0,ps,{qs}),egy,1e-5,"Lennard-Jones Energy");
    compare_vectors(lj.deriv(1,ps,{qs}),grad,1e-5,"Lennard-Jones Gradient");
    compare_vectors(lj.deriv(2,ps,{qs}),hess,1e-5,"Lennard-Jones Hessian");

    test_footer();
    return 0;
} //End main
