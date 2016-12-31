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
#include "ubiquitin.hpp"
#include <cmath>

int main(int argc, char** argv){
    test_header("Testing FourierSeries force-field term");
    FManII::FourierSeries FS;
    
#ifndef NDEBUG
std::vector<double> a(2,0.0),d;
std::map<FManII::Param_t,std::vector<double>> ps={
    {FManII::Param_t::amp,{0.0}},
    {FManII::Param_t::phi,{0.0}},
    {FManII::Param_t::n,{0.0}}};
TEST_THROW(d=FS.deriv(0,ps,{a}),"Qs.size()!=Vs.size()");
ps[FManII::Param_t::amp].push_back(0.0);
TEST_THROW(d=FS.deriv(0,ps,{a}),"Qs.size()!=phis.size()");
ps[FManII::Param_t::phi].push_back(0.0);
TEST_THROW(d=FS.deriv(0,ps,{a}),"Qs.size()!=ns.size()");
#endif
    FManII::CoordArray coords=
        FManII::get_coords(ubiquitin,ubiquitin_conns);
    FManII::ParamSet params=
        FManII::assign_params(coords,FManII::amber99,ubiquitin_FF_types);
    for(auto term:{FManII::IntCoord_t::TORSION,FManII::IntCoord_t::IMPTORSION}){
        const FManII::FFTerm_t term_type=
                std::make_pair(FManII::Model_t::FOURIERSERIES,term);
        const FManII::Vector &qs=coords[term]->get_coords();
        //Energy check
        std::vector<double> Energy=FS.deriv(0,params.at(term_type),{qs});
        if(term==FManII::IntCoord_t::TORSION)
            test_value(Energy[0],ubiquitintorsion_e,1e-5,"Torsion Energy");
        else
            test_value(Energy[0],ubiquitinimproper_e,1e-3,"Improper Torsion Energy");
    }
    
    test_footer();
    return 0;
} //End main


