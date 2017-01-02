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
#include <fstream>
#include "TestMacros.hpp"
#include "testdata/ubiquitin.hpp"
#include "testdata/ubiquitin_params.hpp"

int main(int argc, char** argv){
    test_header("Testing force field parameter assignment function");

    FManII::CoordArray coords=
            FManII::get_coords(ubiquitin,ubiquitin_conns);
    FManII::ParamSet params=
            FManII::assign_params(coords,FManII::amber99,ubiquitin_FF_types);
    for(auto pseti: {std::make_pair(FManII::Model_t::HARMONICOSCILLATOR,
                                    FManII::IntCoord_t::BOND),
                     std::make_pair(FManII::Model_t::HARMONICOSCILLATOR,
                                    FManII::IntCoord_t::ANGLE),
                     std::make_pair(FManII::Model_t::FOURIERSERIES,
                                    FManII::IntCoord_t::TORSION),
                    std::make_pair(FManII::Model_t::FOURIERSERIES,
                                    FManII::IntCoord_t::IMPTORSION),
                    std::make_pair(FManII::Model_t::ELECTROSTATICS,
                                    FManII::IntCoord_t::PAIR14),
                    std::make_pair(FManII::Model_t::ELECTROSTATICS,
                                  FManII::IntCoord_t::PAIR),
                    std::make_pair(FManII::Model_t::LENNARD_JONES,
                                  FManII::IntCoord_t::PAIR14),
                    std::make_pair(FManII::Model_t::LENNARD_JONES,
                                  FManII::IntCoord_t::PAIR),
                     })
    {
        test_value(params.at(pseti),ubiquitin_params.at(pseti),"Params");
    }

    test_footer();
    return 0;

}
