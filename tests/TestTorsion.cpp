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
#include <ForceManII/InternalCoords/Torsion.hpp>
#include "TestMacros.hpp"
#include <cmath>

using namespace std;
using namespace FManII;

int main(int argc, char** argv){
    test_header("Testing torsion internal coordinate");
    vector<Vector> corr_angles({{1.04786541775},{1.04663519795},{3.14110029443},
                                {3.14088699457},{1.0477976969},{1.04666739958},
                                {1.04643698771},{3.14093760341},{1.04778260729}});
    vector<Vector> corr_1st({
        {0.0284653010523,0.97900805365,-0.0302947816383,-0.0464913630787,
         -1.38062723188,-0.185588314476,0.0565934655447,0.916101300855,
         1.0489650447,-0.0385674035183,-0.514482122622,-0.833081948582},
        {-0.0284653010523,-0.97900805365,0.0302947816383,0.033405923143,
         1.36724213856,-0.270576855663,0.00514738517599,-0.852739339491,
         1.10312965726,-0.0100880072668,0.464505254578,-0.86284758324},
        {-0.0284653010523,-0.97900805365,0.0302947816383,0.0284630554296,
         0.97905488329,-0.0304259527828,0.0284818783645,0.978992682927,
         -0.0296820855074,-0.0284796327418,-0.979039512568,0.0298132566519},
        {0.0385758170692,0.515082499748,0.832738175353,-0.0385788138454,
        -0.515253755027,-0.832661578856,-0.0385644067421,-0.514310867344,
         -0.833158545079,0.0385674035183,0.514482122622,0.833081948582},
        {-0.0385758170692,-0.515082499748,-0.832738175353,0.051664253781,
         0.528638848346,1.288826749,-0.0231764439786,0.45094890598,
         -1.31893615688,0.0100880072668,-0.464505254578,0.86284758324},
        {0.0385758170692,0.515082499748,0.832738175353,-0.0566071214945,
         -0.916826103619,-1.04867584612,0.0465109371671,1.38078311644,
         0.186124414111,-0.0284796327418,-0.979039512568,0.0298132566519},
        {-0.0101134825722,0.463945620371,-0.863161242086,0.0232043811408,
         -0.450345967831,1.31929172926,-0.0516583020869,-0.528081775162,
         -1.28921243576,0.0385674035183,0.514482122622,0.833081948582},
        {-0.0101134825722,0.463945620371,-0.863161242086,0.0101189412052,
         -0.46373106115,0.863126559125,0.0100825486338,-0.464719813799,
         0.862882266201,-0.0100880072668,0.464505254578,-0.86284758324},
        {0.0101134825722,-0.463945620371,0.863161242086,-0.00517607349175,
         0.851918316423,-1.10327746201,-0.0334170418223,-1.36701220862,
         0.269929476571,0.0284796327418,0.979039512568,-0.0298132566519}});
    Vector carts({1.1851,-0.0039,0.9875,
                  0.7516,-0.0225,     -0.0209,
                  1.1669,      0.8330,     -0.5693,
                  1.1155,     -0.9329,     -0.5145,
                  -0.7516,      0.0225,      0.0209,
                  -1.1669,    -0.8334,      0.5687,
                  -1.1157,      0.9326,      0.5151,
                  -1.1850,      0.0044,     -0.9875});
    ConnData conns({{1},{0,2,3,4},{1},{1},
                    {1,5,6,7},{4},{4},{4},
                   });
    Torsion d;
    Molecule mol=get_coords(carts,conns);
    auto angles=mol.atom_numbers.at(IntCoord_t::TORSION);
    for(size_t i=0;i<angles.size();++i)
    {
        auto angle=d.deriv(0,carts,angles[i]);
        compare_vectors(angle,corr_angles[i],1e-5,"angles "+to_string(i));

        angle=d.deriv(1,carts,angles[i]);
        compare_vectors(angle,corr_1st[i],1e-5,"1st derivative angle "+to_string(i));

        //bond=d.deriv(2,carts,bonds[i]);
        //compare_vectors(bond,corr_2nd[i],1e-5,"2nd derivative bond "+to_string(i));
    }

    test_footer();
    return 0;
} //End main
