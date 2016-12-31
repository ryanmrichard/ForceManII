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
#include "ubiquitin.hpp"
#include "ubiquitin_intcoords.hpp"


int main(int argc, char** argv){
    test_header("Testing the determination of the internal coordinates");
    
    FManII::CoordArray Coords=FManII::get_coords(ubiquitin,
                                                 ubiquitin_conns);
    for(auto qi:{FManII::IntCoord_t::BOND,
                 FManII::IntCoord_t::ANGLE,
                 FManII::IntCoord_t::TORSION,
                 FManII::IntCoord_t::IMPTORSION,
                 FManII::IntCoord_t::PAIR14,
                 FManII::IntCoord_t::PAIR})
    {
        const std::vector<double>& qs=Coords[qi]->get_coords();
        compare_vectors(qs,ubiquitin_qs.at(qi),1e-4,"Bond displacements");
    }
    test_footer();
    return 0;
}


