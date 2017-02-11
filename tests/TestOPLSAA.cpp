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
#include "testdata/peptide.hpp"
#include "testdata/peptide_deriv.hpp"

using namespace FManII;
using namespace std;

int main(int argc, char** argv){

    test_header("Testing derivatives of pre-packaged OPLS-AA force-fields");
    auto deriv=run_forcemanii(0,peptide,peptide_conns,oplsaa,
                                    peptide_FF_types);
    for(auto derivi : deriv)
    {
        const auto ffterm=derivi.first;
        const string msg=string("OPLS-AA")+" "+ffterm.first+" "+ffterm.second+" energy";
        double tol=1e-5;
        if(ffterm.first==Model_t::ELECTROSTATICS||
           ffterm.second==IntCoord_t::IMPTORSION)tol=6e-4;
        if(ffterm.second==IntCoord_t::PAIR13)continue;//It's zero
        test_value(derivi.second[0],peptide_egys.at(ffterm),tol,msg);
    }

    deriv=run_forcemanii(1,peptide,peptide_conns,oplsaa,
                                    peptide_FF_types);
    for(auto derivi : deriv)
    {
        const auto ffterm=derivi.first;
        if(ffterm.second==IntCoord_t::PAIR13)continue;//It's zero
        const string msg=string("OPLS-AA")+" "+ffterm.first+" "+ffterm.second+" gradient";
        compare_vectors(derivi.second,corr_grad.at(ffterm),1e-5,msg);
    }

    test_footer();
    return 0;
}
