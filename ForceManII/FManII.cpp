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

namespace FManII {

CoordArray get_coords(const Vector& Carts,
                      const ConnData& Conns){
    const size_t NAtoms=Carts.size()/3;
    DEBUG_CHECK(NAtoms==Conns.size(),"Number of atoms differs among inputs");
    auto Sys=std::make_shared<Vector>(Carts);
    
    CoordArray FoundCoords;
    FoundCoords.emplace(IntCoord_t::BOND,std::move(make_unique<Bond>(Sys)));
    FoundCoords.emplace(IntCoord_t::ANGLE,std::move(make_unique<Angle>(Sys)));
    FoundCoords.emplace(IntCoord_t::TORSION,std::move(make_unique<Torsion>(Sys)));
    FoundCoords.emplace(IntCoord_t::IMPTORSION,std::move(make_unique<ImproperTorsion>(Sys)));
    FoundCoords.emplace(IntCoord_t::PAIR,std::move(make_unique<Pair>(Sys)));
    FoundCoords.emplace(IntCoord_t::PAIR13,std::move(make_unique<Pair13>(Sys)));
    FoundCoords.emplace(IntCoord_t::PAIR14,std::move(make_unique<Pair14>(Sys)));
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
                    FoundCoords[IntCoord_t::TORSION]->add_coord({AtomI,AtomJ,AtomK,AtomL});
                }//Close AtomL
                if(AtomK<AtomI)continue;//Avoid 2x counting angle
                pair13.insert(std::make_pair(AtomI,AtomK));
                FoundCoords[IntCoord_t::ANGLE]->add_coord({AtomI,AtomJ,AtomK});
                if(Conns[AtomJ].size()==3){
                    for(size_t AtomL: Conns[AtomJ]){
                        if(AtomL==AtomK||AtomL==AtomI||AtomL<AtomK)continue;
                        FoundCoords[IntCoord_t::IMPTORSION]->add_coord({AtomI,AtomJ,AtomK,AtomL});
                    }//Close AtomL imptorsion
                }//Close if imptorsion
            }//Close Atom K
            if(AtomJ<AtomI)continue; //Avoid 2x counting bond
            pair12.insert(std::make_pair(AtomI,AtomJ));
            FoundCoords[IntCoord_t::BOND]->add_coord({AtomI,AtomJ});
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
            FoundCoords[ctype]->add_coord({AtomI,AtomJ});
        }
    }
    return FoundCoords;
}

ParamSet assign_params(const CoordArray& coords,
                       const ForceField& ff,
                       const IVector& Types,
                       bool skip_missing)
 {
    ParamSet ps;
    for(const auto& termi:ff.terms){
        const FFTerm_t term_type=termi.first;
        for(const auto& coordi:termi.second.coords){
            const InternalCoordinates& coord=*coords.at(coordi);
            for(auto parami:termi.second.model().params){
                ps[term_type][parami]=
                    ff.assign_param(term_type,parami,coord,Types,skip_missing);
            }
        }
    }
    return ps;
}

DerivType deriv(size_t order,
                                const ForceField& ff,
                                const ParamSet& ps,
                                const CoordArray& coords)
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

HarmonicBond::HarmonicBond():
    FFTerm(std::make_shared<HarmonicOscillator>(),{IntCoord_t::BOND}){

}

HarmonicAngle::HarmonicAngle():
    FFTerm(std::make_shared<HarmonicOscillator>(),{IntCoord_t::ANGLE}){}

FourierTorsion::FourierTorsion():
    FFTerm(std::make_shared<FourierSeries>(),{IntCoord_t::TORSION}){}

FourierImproperTorsion::FourierImproperTorsion():
    FFTerm(std::make_shared<FourierSeries>(),{IntCoord_t::IMPTORSION}){}

LJ14::LJ14():
    FFTerm(std::make_shared<LennardJones>(),{IntCoord_t::PAIR14}){}

LJPair::LJPair():
    FFTerm(std::make_shared<LennardJones>(),{IntCoord_t::PAIR}){}

Electrostatics14::Electrostatics14():
    FFTerm(std::make_shared<Electrostatics>(),{IntCoord_t::PAIR14}){}

ElectrostaticsPair::ElectrostaticsPair():
    FFTerm(std::make_shared<Electrostatics>(),{IntCoord_t::PAIR}){}

} //End namespace FManII

