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

#include "FManIIPulsarAPI.hpp"
#include "ForceManII/FManII.hpp"
#include <iostream>
using namespace FManII;
using namespace std;
using namespace pulsar;


DerivReturnType FFTermPulsar::deriv_(size_t Order,const Wavefunction& wfn)
{
    auto ff_name=options().get<string>("FORCE_FIELD");
    auto types=options().get<vector<size_t>>("ATOM_TYPES");
    auto model_name=options().get<string>("MODEL_NAME");
    auto intcoord_name=options().get<string>("COORD_NAME");    Conn_t temp_conns=get_connectivity(*wfn.system);
    vector<Atom> index2atom;
    unordered_map<Atom,size_t> atom2index;
    for(auto atomi:*wfn.system)
    {
        atom2index[atomi]=index2atom.size();
        index2atom.push_back(atomi);
    }
    ConnData conns(index2atom.size());
    Vector Carts;
    for(size_t i=0;i<index2atom.size();++i)
    {
        for(size_t j=0;j<3;j++)
            Carts.push_back(index2atom[i][j]);
        for(auto atomj:temp_conns.at(index2atom[i]))
            conns[i].insert(atom2index.at(atomj));
    }
    auto deriv_comps=run_forcemanii(Order,Carts,conns,get_ff(ff_name),types);
    if(model_name=="ALL" && intcoord_name=="ALL")
    {
       Vector deriv(std::pow(Carts.size(),Order));
       for(const auto& dci:deriv_comps){
           const auto& di=dci.second;
           for(size_t i=0;i<di.size();++i)
               deriv[i]+=di[i];
       }
       return {wfn,deriv};
    }
    return {wfn,deriv_comps.at(make_pair(model_name,intcoord_name))};
}
