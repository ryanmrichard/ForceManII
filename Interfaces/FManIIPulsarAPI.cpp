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
#include "ForceManII/Common.hpp"
#include <iostream>
using namespace FManII;
using namespace std;
using namespace pulsar;


pair<Vector,ConnData> make_carts_conns(const Wavefunction& wfn)
{
    Conn_t temp_conns=get_connectivity(*wfn.system);
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
    return {Carts,conns};
}

DerivReturnType FFPulsar::deriv_(size_t Order,const Wavefunction& wfn)
{
    auto ff_name=options().get<string>("FORCE_FIELD");
    auto types=options().get<vector<size_t>>("ATOM_TYPES");
    auto m2s=options().get<vector<string>>("SKIP_MODELS");
    auto c2s=options().get<vector<string>>("SKIP_COORDS");
    CHECK(types.size()==wfn.system->size(),
          "Number of atom types must match number of atoms "+
          to_string(types.size())+":"+to_string(wfn.system->size()));
    pair<Vector,ConnData> geom=make_carts_conns(wfn);
    const Vector& Carts=geom.first;
    const ConnData& conns=geom.second;
    auto deriv_comps=run_forcemanii(Order,Carts,conns,get_ff(ff_name),types);
    Vector deriv(std::pow(Carts.size(),Order));
    for(const auto& dci:deriv_comps){
        if(find(m2s.begin(),m2s.end(),dci.first.first)!=m2s.end() ||
           find(c2s.begin(),c2s.end(),dci.first.second)!=c2s.end())continue;
        const auto& di=dci.second;
        for(size_t i=0;i<di.size();++i)
            deriv[i]+=di[i];
    }
    return {wfn,deriv};
}


DerivReturnType FFTermPulsar::deriv_(size_t Order,const Wavefunction& wfn)
{
    auto ff_name=options().get<string>("FORCE_FIELD");
    auto types=options().get<vector<size_t>>("ATOM_TYPES");
    auto model_name=options().get<string>("MODEL_NAME");
    auto intcoord_name=options().get<string>("COORD_NAME");
    CHECK(types.size()==wfn.system->size(),
          "Number of atom types must match number of atoms "+
          to_string(types.size())+":"+to_string(wfn.system->size()));
    pair<Vector,ConnData> geom=make_carts_conns(wfn);
    const Vector& Carts=geom.first;
    const ConnData& conns=geom.second;
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

Vector FFCharges::calculate_(unsigned int deriv,
                             const Wavefunction & wfn){
    auto ff_name=options().get<string>("FORCE_FIELD");
    auto types=options().get<vector<size_t>>("ATOM_TYPES");
    const size_t natoms=types.size();
    const size_t DoF=std::pow(3*natoms,deriv);
    if(deriv>0)
        return Vector(0.0,natoms*DoF);

    pair<Vector,ConnData> geom=make_carts_conns(wfn);
    const Vector& Carts=geom.first;
    const ConnData& conns=geom.second;
    auto ff=get_ff(ff_name);
    Vector qs;
    for(size_t i : types)
        qs.push_back(ff.params.get_param(Terms_t::CL,Param_t::q,{i})[0]);
    return qs;

}
