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

using namespace std;
using namespace FManII;
const map<string,ForceField> name2ff({{"amber99",amber99},
                                      {"oplsaa",oplsaa},
                                      {"charmm22",charmm22}});

int main(int argc, char** argv){
    test_header("Testing force field parsing function");
    for(auto ffs:name2ff){
        const string &ffname=ffs.first;
        ifstream file("../../ForceFields/"+ffname+".prm");
        ForceField ff=parse_file(move(file),ffs.first=="charmm22");
        ForceField corr_ff=ffs.second;
        for(auto pi: corr_ff.params){//Compare each parameter in the FF
            const string &fname=get<0>(pi).first,//FFTerm
                         &cname=get<0>(pi).second,//IntCoord
                         &pname=get<1>(pi);//Parameter type
            const IVector &types=get<2>(pi);
            string stypes="";
            for(size_t s : types)stypes+=to_string(s)+" ";
            string msg=ffname+" "+fname+" "+cname+" "+stypes+" "+pname;
            const auto& da_param=ff.params.get_param(get<0>(pi),pname,types);
            compare_vectors(da_param,get<3>(pi),1e-10,msg);
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
