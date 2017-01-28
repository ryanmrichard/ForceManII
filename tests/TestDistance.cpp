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
#include <ForceManII/InternalCoords/Distance.hpp>
#include "TestMacros.hpp"
#include <cmath>

using namespace std;
using namespace FManII;

int main(int argc, char** argv){
    test_header("Testing distance internal coordinate");
    vector<Vector> corr_distance({{0.855862138431},{0.74672619}});
    vector<Vector> corr_1st({
      {-0.99315060432287627, -0.11684124756739722, 0.0,
       0.99315060432287627, 0.11684124756739722, -0.0},
      {-0.1339179,-0.99099243,0.,
        0.1339179, 0.99099243, -0.}
    });

    Vector carts({0.0,0.0,0.0,0.85,0.1,0.0,0.1,0.74,0.0});
    ConnData conns({{1,2},{0},{0}});
    IVector types({2001,2002,2002});
    Bond d;
    Molecule mol=get_coords(carts,conns);
    auto bonds=mol.atom_numbers.at(IntCoord_t::BOND);
    for(size_t i=0;i<bonds.size();++i)
    {
        auto bond=d.deriv(0,carts,bonds[i]);
        compare_vectors(bond,corr_distance[i],1e-5,"distance "+to_string(i));

        bond=d.deriv(1,carts,bonds[i]);
        compare_vectors(bond,corr_1st[i],1e-5,"1st derivative bond "+to_string(i));
    }

    test_footer();
    return 0;
} //End main
