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
#include "testdata/ubiquitin_intcoords.hpp"
#include "testdata/peptide.hpp"
#include "testdata/peptide_intcoords.hpp"

using namespace std;
using namespace FManII;

int main(int argc, char** argv){
    test_header("Testing the determination of the internal coordinates");

    for(auto mol:{make_tuple(ubiquitin,ubiquitin_conns,ubiquitin_qs,"ubiquitin"),
                  make_tuple(peptide,peptide_conns,peptide_qs,"peptide")})
    {
        CoordArray Coords=get_coords(get<0>(mol),get<1>(mol));
        for(auto qi:{IntCoord_t::BOND,IntCoord_t::ANGLE,
                     IntCoord_t::TORSION,IntCoord_t::IMPTORSION,
                     IntCoord_t::PAIR14,IntCoord_t::PAIR})
        {
            const auto& qs=Coords[qi]->get_coords();
            compare_vectors(qs,get<2>(mol).at(qi),1e-4,string(get<3>(mol))+" "+qi);
        }
    }
    test_footer();
    return 0;
}


