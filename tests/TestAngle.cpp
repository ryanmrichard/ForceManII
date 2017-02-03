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
#include <ForceManII/InternalCoords/Angle.hpp>
#include "TestMacros.hpp"
#include <cmath>

using namespace std;
using namespace FManII;

int main(int argc, char** argv){
    test_header("Testing angle internal coordinate");
    vector<Vector> corr_angles({{1.31936614028}});
    vector<Vector> corr_1st({
      {0.13651877133105803,-1.1604095563139933, 0.0,
       1.190597441007536, 0.98106952761958865, 0.0,
       -1.3271162123385942, 0.1793400286944046, 0.0}
    });
    vector<Vector> corr_2nd({
       {},
       {}
    });

    Vector carts({0.0,0.0,0.0,0.85,0.1,0.0,0.1,0.74,0.0});
    ConnData conns({{1,2},{0},{0}});
    IVector types({2001,2002,2002});
    Angle d;
    Molecule mol=get_coords(carts,conns);
    auto angles=mol.atom_numbers.at(IntCoord_t::ANGLE);
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
