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

int main(int argc, char** argv){
    test_header("Testing Distance class");
 
    FManII::CoordArray Bonds=FManII::get_coords(ubiquitin,ubiquitin_FF_types,
         ubiquitin_FF_params,ubiquitin_conns,1/1.2);
    const std::vector<double>& BondLength=
        Bonds[FManII::IntCoord_t::BOND]->values();
    const std::vector<double>& bond_k=Bonds[FManII::BOND]->params(FManII::K);
    
    compare_vectors(BondLength,ubiquitin_bonds,1e-4,"Bond displacements");
    compare_vectors(bond_k,ubiquitin_K,1e-4,"Bond parameters");
    
    const std::vector<double>& PairLength=
        Bonds[FManII::ELECTROSTATICS]->values();
    const std::vector<double>& Pair_q=
        Bonds[FManII::ELECTROSTATICS]->params(FManII::q);
    compare_vectors(PairLength,ubiquitin_pairs,1e-4,"Pair distances");
    compare_vectors(Pair_q,ubiquitin_charge_q,1e-4,"Charges");
    
    
    test_footer();
    return 0;
}

