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
#include "TestMacros.hpp"
#include "testdata/ubiquitin.hpp"
#include "testdata/peptide.hpp"
#include <cmath>

using namespace FManII;
using namespace std;

int main(int argc, char** argv){

    test_header("Testing derivatives of pre-packaged force-fields");
    auto imp=make_pair<string,string>(Model_t::FOURIERSERIES,
                                      IntCoord_t::IMPTORSION);

    auto ub_amber=make_pair(run_forcemanii(0,ubiquitin,ubiquitin_conns,amber99,
                                           ubiquitin_FF_types),"AMBER99");
    auto pep_oplsaa=make_pair(run_forcemanii(0,peptide,peptide_conns,oplsaa,
                                           peptide_FF_types), "OPLSAA");
    for(auto ff:{ub_amber,pep_oplsaa})
    {
        for(auto derivi : ff.first)
        {
            const auto ffterm=derivi.first;
            const string msg=string(ff.second)+" "+ffterm.first+" "+
                                  ffterm.second+" energy";
            double tol=1e-5;
            if(ffterm.first==Model_t::ELECTROSTATICS||ffterm==imp)tol=4e-4;
            if(ff.second=="AMBER99")
                test_value(derivi.second[0],ubiquitin_egys.at(ffterm),tol,msg);
            else if(ff.second=="OPLSAA")
                test_value(derivi.second[0],peptide_egys.at(ffterm),tol,msg);
        }
    }

    test_footer();
    return 0;
}
