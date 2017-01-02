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

int main(int argc, char** argv){
    test_header("Testing FourierSeries force-field term");
    FManII::FourierSeries FS;
    
#ifndef NDEBUG
std::vector<double> a(2,0.0),d;
std::map<std::string,std::vector<double>> ps={
    {FManII::Param_t::amp,{0.0}},
    {FManII::Param_t::phi,{0.0}},
    {FManII::Param_t::n,{0.0}}};
TEST_THROW(d=FS.deriv(0,ps,{a}),"Qs.size()!=Vs.size()");
ps[FManII::Param_t::amp].push_back(0.0);
TEST_THROW(d=FS.deriv(0,ps,{a}),"Qs.size()!=phis.size()");
ps[FManII::Param_t::phi].push_back(0.0);
TEST_THROW(d=FS.deriv(0,ps,{a}),"Qs.size()!=ns.size()");
#endif
    FManII::DerivType deriv=FManII::run_forcemanii(0,ubiquitin,ubiquitin_conns,
                                           FManII::amber99,ubiquitin_FF_types);

    for(auto term:{FManII::IntCoord_t::TORSION,FManII::IntCoord_t::IMPTORSION}){
        const FManII::FFTerm_t term_type=
                std::make_pair(FManII::Model_t::FOURIERSERIES,term);
        const double tol=term==FManII::IntCoord_t::TORSION?1e-5:3e-4;
        test_value(deriv.at(term_type)[0],ubiquitin_egys.at(term_type),tol,"Torsion Energy");
    }
    
    test_footer();
    return 0;
} //End main


