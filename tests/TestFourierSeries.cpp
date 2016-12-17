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
#include <ForceManII/FourierSeries.hpp>
#include "TestMacros.hpp"
#include "ubiquitin.hpp"
#include <cmath>

int main(int argc, char** argv){
    test_header("Testing FourierSeries force-field term");
    FManII::FourierSeries FS;
    
#ifndef NDEBUG
std::vector<double> a(3,0.0),b(2,0.0),c,d;
TEST_THROW(d=FS.deriv(a,c,b),"Qs.size()!=Vs.size()");
TEST_THROW(d=FS.deriv(a,b,c),"Qs.size()!=ns.size()");
#endif

     FManII::CoordArray coords=FManII::get_coords(ubiquitin,ubiquitin_FF_types,
            ubiquitin_FF_params,ubiquitin_conns);
     const std::vector<double>& torsions=coords[FManII::TORSION]->values(),
        tor_v=coords[FManII::TORSION]->params(FManII::amp),
        tor_n=coords[FManII::TORSION]->params(FManII::n),
        imps=coords[FManII::IMPTORSION]->values(),
        imp_v=coords[FManII::IMPTORSION]->params(FManII::amp),
        imp_n=coords[FManII::IMPTORSION]->params(FManII::n);
     
    //Energy check
    std::vector<double> Energy=FS.deriv(torsions,tor_v,tor_n);
    test_value(Energy[0],ubiquitintorsion_e,1e-5,"Torsion Energy");
    
    Energy=FS.deriv(imps,imp_v,imp_n);
    
    //Lower convergence than other terms b/c of ambiguity in imp
    test_value(Energy[0],ubiquitinimproper_e,1e-3,"Improper Torsion Energy");
    
    test_footer();
    return 0;
} //End main


