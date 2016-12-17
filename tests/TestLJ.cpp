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
#include <ForceManII/LennardJones.hpp>
#include "TestMacros.hpp"
#include "ubiquitin.hpp"
#include <cmath>

int main(int argc, char** argv){
    test_header("Testing 6-12 force-field term");
    FManII::LennardJones lj;

#ifndef NDEBUG
std::vector<double> a(3,0.0),b(2,0.0),c;
TEST_THROW(c=lj.deriv(a,b,a),"sigmas.size()==dist.size()");
TEST_THROW(c=lj.deriv(a,a,b),"epsilons.size()==dist.size()");
#endif

     FManII::CoordArray coords=FManII::get_coords(ubiquitin,ubiquitin_FF_types,
            ubiquitin_FF_params,ubiquitin_conns,1/1.2,1/2.0);
     const std::vector<double>& dist=
        coords[FManII::LENNARD_JONES]->values();
     const std::vector<double>& sigmas=
        coords[FManII::LENNARD_JONES]->params(FManII::sigma);
     const std::vector<double>& epsilons=
        coords[FManII::LENNARD_JONES]->params(FManII::epsilon);

    //Energy check
    std::vector<double> Energy=lj.deriv(dist,sigmas,epsilons);
    test_value(Energy[0],ubiquitinvdw_e,1e-5,"Lennard-Jones 6-12 Energy");

    test_footer();
    return 0;
} //End main
