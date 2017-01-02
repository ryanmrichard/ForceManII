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
    for(std::string ffname:{"amber99"}){
        std::ifstream file("../../ForceFields/"+ffname+".prm");
        FManII::ForceField ff=FManII::parse_file(std::move(file));
        FManII::ForceField corr_ff=
                ffname=="amber99"?FManII::amber99:FManII::amber99;
        //Loops over sets of parameters (e.g. harmonic bond-stretching Ks) for
        //finer error highlighting
        for(const auto& model:ff.params)
            for(const auto& param:model.second){
                const auto& actual=
                    corr_ff.params.at(model.first).at(param.first);
            test_value(param.second,actual,ffname+" params");
        }
        test_value(ff.type2class,corr_ff.type2class,ffname+" Type2Class");
        test_value(ff.terms,corr_ff.terms,ffname+" terms");
        test_value(ff.orderrules,corr_ff.orderrules,ffname+" order rules");
        test_value(ff.paramtypes,corr_ff.paramtypes,ffname+"p types");
        test_value(ff.combrules,corr_ff.combrules,ffname+" combrules");
        test_value(ff.scale_factors,corr_ff.scale_factors,ffname+" scale");
        test_value(ff,corr_ff,ffname+" parsed successfully");
    }

    test_footer();
    return 0;
}
