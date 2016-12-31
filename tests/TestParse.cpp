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

int main(int argc, char** argv){
    test_header("Testing force field parsing function");

    std::ifstream file("../../ForceFields/amber99.prm");
    FManII::ForceField ff=FManII::parse_file(std::move(file));
    //Loops over sets of parameters (e.g. harmonic bond-stretching Ks) for
    //finer error highlighting
    for(const auto& model:ff.params)
        for(const auto& param:model.second){
            const auto& actual=
                    FManII::amber99.params.at(model.first).at(param.first);
            test_value(param.second,actual,"Amber99 params");
        }
    test_value(ff.type2class,FManII::amber99.type2class,"Amber99 Type2Class");
    test_value(ff.terms,FManII::amber99.terms,"Amber99 terms");
    test_value(ff.orderrules,FManII::amber99.orderrules,"Amber order rules");
    test_value(ff.paramtypes,FManII::amber99.paramtypes,"Amber99 p types");
    test_value(ff.combrules,FManII::amber99.combrules,"Amber99 combrules");
    test_value(ff.scale_factors,FManII::amber99.scale_factors,"Amber99 scale");
    test_value(ff,FManII::amber99,"Amber99.prm parsed successfully");

    test_footer();
    return 0;
}
