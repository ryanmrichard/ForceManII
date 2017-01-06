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

using namespace std;
using namespace FManII;

int main(int argc, char** argv){
    test_header("Testing charge-charge force-field term");
    Electrostatics ct;
    Vector Qs({3.0,2.0,1.0}),c;
    map<string,Vector> qs({{Param_t::q,{0.0}}});
#ifndef NDEBUG
TEST_THROW(c=ct.deriv(0,qs,{Qs}),"charges.size()==dist.size()");
TEST_THROW(c=ct.deriv(4,qs,{Qs}),"Derivative higher than 2 requested");
#endif
    qs[Param_t::q]=Vector({6.0,7.0,8.0});
    const Vector egy={13.5},grad={-0.66666667,-1.75,-8.},
          hess={0.44444444,0.0,0.0,0.0,1.75,0.0,0.0,0.0,16.0};
    compare_vectors(ct.deriv(0,qs,{Qs}),egy,1e-5,"Coulomb's law energy");
    compare_vectors(ct.deriv(1,qs,{Qs}),grad,1e-5,"Gradient of Coulomb's law");
    compare_vectors(ct.deriv(2,qs,{Qs}),hess,1e-5,"Hessian of Coulomb's law");
        
    test_footer();
    return 0;
} //End main
