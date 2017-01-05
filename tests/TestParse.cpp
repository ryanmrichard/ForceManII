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
    for(std::string ffname:{"amber99","oplsaa"}){
        std::ifstream file("../../ForceFields/"+ffname+".prm");
        FManII::ForceField ff=FManII::parse_file(std::move(file));
        FManII::ForceField corr_ff=
                ffname=="amber99"?FManII::amber99:FManII::oplsaa;
        for(auto pi: ff.params){//Compare each parameter in the FF
            const std::string fname=std::get<0>(pi).first,//FFTerm
                              cname=std::get<0>(pi).second;//IntCoord
            const std::string pname=std::get<1>(pi);//Parameter type
            FManII::IVector types=std::get<2>(pi);
            std::string stypes="";
            for(size_t s : types)stypes+=std::to_string(s)+" ";
            std::string msg=ffname+" "+fname+" "+cname+" "+stypes+" "+pname;
            const auto& corr_param=
                    corr_ff.params.get_param(std::get<0>(pi),pname,types);
            compare_vectors(std::get<3>(pi),corr_param,1e-10,msg);
        }
        test_value(ff.type2class,corr_ff.type2class,ffname+" Type2Class");
        test_value(ff.terms,corr_ff.terms,ffname+" terms");
        test_value(ff.orderrules,corr_ff.orderrules,ffname+" order rules");
        test_value(ff.paramtypes,corr_ff.paramtypes,ffname+"p types");
        test_value(ff.combrules,corr_ff.combrules,ffname+" combrules");
        test_value(ff.scale_factors,corr_ff.scale_factors,ffname+" scale");
    }

    test_footer();
    return 0;
}
