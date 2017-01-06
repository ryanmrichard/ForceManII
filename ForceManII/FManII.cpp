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
#include "ForceManII/FManII.hpp"
#include "ForceManII/Common.hpp"
#include "ForceManII/Util.hpp"
#include "ForceManII/InternalCoords/Distance.hpp"
#include "ForceManII/InternalCoords/Angle.hpp"
#include "ForceManII/InternalCoords/Torsion.hpp"
#include "ForceManII/InternalCoords/ImproperTorsion.hpp"
#include <iostream>
#include <algorithm>

using namespace std;
namespace FManII {

inline void add_coord(Molecule& FoundCoords,
                      const std::string& name,
                      const IVector& atoms)
{
    FoundCoords.atom_numbers[name].push_back(atoms);
    FoundCoords.coords[name].push_back(
            get_intcoord(name)->deriv(0,*FoundCoords.carts,atoms)[0]);

}


Molecule get_coords(const Vector& Carts,
                      const ConnData& Conns){
    const size_t NAtoms=Carts.size()/3;
    DEBUG_CHECK(NAtoms==Conns.size(),"Number of atoms differs among inputs");
    auto Sys=std::make_shared<Vector>(Carts);
    
    Molecule FoundCoords;
    FoundCoords.carts=Sys;
    std::set<std::pair<size_t,size_t>> pair14,pair13,pair12;
    for(size_t AtomI=0;AtomI<NAtoms;++AtomI){
        for(size_t AtomJ : Conns[AtomI]){
            DEBUG_CHECK(AtomI!=AtomJ,"AtomI is bonded to itself. What??");
            for(size_t AtomK: Conns[AtomJ]){
                DEBUG_CHECK(AtomK!=AtomJ,"Atom J is bonded to itself. What??");
                if(AtomK==AtomI)continue;//Went backwards in the graph
                for(size_t AtomL : Conns[AtomK]){
                    DEBUG_CHECK(AtomL!=AtomK,"Atom K is bonded to itself");
                    if(AtomL==AtomJ)continue;//Went backwards
                    if(AtomL>AtomI){
                        pair14.insert(std::make_pair(AtomI,AtomL));
                    }
                    if(AtomK<AtomJ)continue;
                    add_coord(FoundCoords,IntCoord_t::TORSION,{AtomI,AtomJ,AtomK,AtomL});
                }//Close AtomL
                if(AtomK<AtomI)continue;//Avoid 2x counting angle
                if(!Conns[AtomI].count(AtomK))//Ensure it's not a three-membered ring
                {
                    pair13.insert(std::make_pair(AtomI,AtomK));
                    add_coord(FoundCoords,IntCoord_t::PAIR13,{AtomI,AtomJ,AtomK});
                }
                add_coord(FoundCoords,IntCoord_t::ANGLE,{AtomI,AtomJ,AtomK});
                if(Conns[AtomJ].size()==3){
                    for(size_t AtomL: Conns[AtomJ]){
                        if(AtomL==AtomK||AtomL==AtomI||AtomL<AtomK)continue;
                        add_coord(FoundCoords,IntCoord_t::IMPTORSION,{AtomI,AtomJ,AtomK,AtomL});
                    }//Close AtomL imptorsion
                }//Close if imptorsion
            }//Close Atom K
            if(AtomJ<AtomI)continue; //Avoid 2x counting bond
            pair12.insert(std::make_pair(AtomI,AtomJ));
            add_coord(FoundCoords,IntCoord_t::BOND,{AtomI,AtomJ});
        }//Close Atom J     
    }//Close AtomI
    for(size_t AtomI=0;AtomI<NAtoms;++AtomI){
        for(size_t AtomJ=AtomI+1;AtomJ<NAtoms;++AtomJ){
            auto pair=std::make_pair(AtomI,AtomJ);
            if(pair12.count(pair))continue;
            const bool is13=pair13.count(pair),is14=pair14.count(pair);
            std::string ctype=IntCoord_t::PAIR;
            if(is13)ctype=IntCoord_t::PAIR13;
            else if(is14)ctype=IntCoord_t::PAIR14;
            add_coord(FoundCoords,ctype,{AtomI,AtomJ});
        }
    }
    return FoundCoords;
}

ParamSet assign_params(const Molecule& sys,
                       const ForceField& ff,
                       const IVector& Types,
                       bool skip_missing)
 {
    ParamSet ps;
    for(const auto& termi:ff.terms){
        const FFTerm_t term_type=termi.first;
        //for(const auto& coordi:termi.second.coords){
            for(auto parami:termi.second.model().params){
                ps[term_type][parami]=
                    ff.assign_param(term_type,parami,sys.atom_numbers.at(term_type.second),Types,skip_missing);
            }
        //}
    }
    return ps;
}

DerivType deriv(size_t order,
                const ForceField& ff,
                const ParamSet& ps,
                const Molecule& coords)
{
   DerivType rv;
   for(const auto& i:ff.terms){
        const auto& term_type=i.first;
        const FFTerm& ffterm=i.second;
        Vector d=ffterm.deriv(order,ps.at(term_type),coords);
        if(ff.scale_factors.count(term_type))
            for(double& di:d)di*=ff.scale_factors.at(term_type);
        rv.emplace(term_type,std::move(d));
    }
    return rv;
}

shared_ptr<ModelPotential> get_potential(const string& name)
{
    if(name==Model_t::HARMONICOSCILLATOR)
        return make_shared<HarmonicOscillator>();
    else if(name==Model_t::FOURIERSERIES)
        return make_shared<FourierSeries>();
    else if(name==Model_t::ELECTROSTATICS)
        return make_shared<Electrostatics>();
    else if(name==Model_t::LENNARD_JONES)
        return make_shared<LennardJones>();
    throw runtime_error(name+" is not a known model potential");
}

shared_ptr<InternalCoordinates> get_intcoord(const string& name)
{
    if(name==IntCoord_t::BOND)
        return make_shared<Bond>();
    else if(name==IntCoord_t::ANGLE)
        return make_shared<Angle>();
    else if(name==IntCoord_t::TORSION)
        return make_shared<Torsion>();
    else if(name==IntCoord_t::IMPTORSION)
        return make_shared<ImproperTorsion>();
    else if (name==IntCoord_t::PAIR)
        return make_shared<Pair>();
    else if(name==IntCoord_t::PAIR13)
        return make_shared<Pair13>();
    else if(name==IntCoord_t::PAIR14)
        return make_shared<Pair14>();
    throw runtime_error(name+" is not a known internal coordinate");
}

} //End namespace FManII

