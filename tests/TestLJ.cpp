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
#include "ubiquitin.hpp"
#include <cmath>

int main(int argc, char** argv){
    test_header("Testing 6-12 force-field term");
    FManII::LennardJones lj;

#ifndef NDEBUG
std::map<FManII::Param_t,std::vector<double>> ps={
    {FManII::Param_t::sigma,{0.0}},
    {FManII::Param_t::epsilon,{0.0}}};
std::vector<double> a(2,0.0),c;
TEST_THROW(c=lj.deriv(0,ps,{a}),"Wrong number of parameters");
TEST_THROW(c=lj.deriv(0,ps,{a,a}),"Wrong number of coordinates");
TEST_THROW(c=lj.deriv(0,ps,{a}),"sigmas.size()==dist.size()");
ps[FManII::Param_t::sigma].push_back(0.0);
TEST_THROW(c=lj.deriv(0,ps,{a}),"epsilons.size()==dist.size()");
#endif

    FManII::CoordArray coords=
        FManII::get_coords(ubiquitin,ubiquitin_conns);
    FManII::ParamSet params=
        FManII::assign_params(coords,FManII::amber99,ubiquitin_FF_types);
    std::vector<double> Energy(1,0.0);
    for(auto param:{FManII::IntCoord_t::PAIR,FManII::IntCoord_t::PAIR14}){
        const std::vector<double> &dist=coords[param]->get_coords();
        auto term_type=std::make_pair(FManII::Model_t::ELECTROSTATICS,param);
        double scale=1.0;
        if(FManII::amber99.scale_factors.count(term_type))
            scale=FManII::amber99.scale_factors.at(term_type);
        Energy[0]+=scale*lj.deriv(0,params[term_type],{dist})[0];
    }

    test_value(Energy[0],ubiquitinvdw_e,1e-5,"Lennard-Jones 6-12 Energy");

    test_footer();
    return 0;
} //End main
