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
#include "ForceManII/Distance.hpp"
#include "ForceManII/Angle.hpp"
#include "ForceManII/Torsion.hpp"
#include "ForceManII/Common.hpp"
#include "ForceManII/Util.hpp"
#include <iostream>
#include <algorithm>

namespace FManII {
using DVector=std::vector<double>;
using IVector=std::vector<size_t>;
template<size_t n> using IArray=std::array<size_t,n>;
using shared_DVector=std::shared_ptr<const DVector>;
using std::max;
using std::min;

/* These end up being giant messy blocks of indirection that disguises what's
   actually going on so we introduce some functions to hide it and make it
   prettier to read*/
inline void add_bond(const IVector& Atoms,const AtomTypes& Types,
                     const ParamTypes& Params,CoordArray& FoundCoords);

inline void add_angle(const IVector& Atoms,
               const AtomTypes& Types,const ParamTypes& Params,
               CoordArray& FoundCoords);

inline void add_torsion(const IVector& Atoms,
        const AtomTypes& Types,const ParamTypes& Params,
        CoordArray& FoundCoords);

inline void add_imp(const IVector& Atoms,
        const AtomTypes& Types,const ParamTypes& Params,
        CoordArray& FoundCoords,const DVector& carts);

inline void add_els(const IVector& Atoms,const AtomTypes& Types,
        const ParamTypes& Params,CoordArray& FoundCoords,double scale);
inline void add_lj(const IVector& Atoms,const AtomTypes& Types,
        const ParamTypes& Params,CoordArray& FoundCoords,double scale);

CoordArray get_coords(const DVector& Carts,const AtomTypes& Types,
                     const ParamTypes& Params,const ConnData& Conns,
                     double chg14scale,double vdw14scale){
    const size_t NAtoms=Carts.size()/3;
    DEBUG_CHECK(NAtoms==Conns.size(),"Number of atoms differs among inputs");
    shared_DVector Sys=std::make_shared<DVector>(Carts);
    
    CoordArray FoundCoords;
    FoundCoords.emplace(BOND,std::move(make_unique<Distance>(Sys)));
    FoundCoords.emplace(ANGLE,std::move(make_unique<Angle>(Sys)));
    FoundCoords.emplace(TORSION,std::move(make_unique<Torsion>(Sys)));
    FoundCoords.emplace(IMPTORSION,std::move(make_unique<Torsion>(Sys)));
    FoundCoords.emplace(ELECTROSTATICS,std::move(make_unique<Distance>(Sys)));
    FoundCoords.emplace(LENNARD_JONES,std::move(make_unique<Distance>(Sys)));
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
                    add_torsion({AtomI,AtomJ,AtomK,AtomL},Types,Params,FoundCoords);        
                }//Close AtomL
                if(AtomK<AtomI)continue;//Avoid 2x counting angle
                pair13.insert(std::make_pair(AtomI,AtomK));
                add_angle({AtomI,AtomJ,AtomK},Types,Params,FoundCoords);
                if(Conns[AtomJ].size()==3){
                    for(size_t AtomL: Conns[AtomJ]){
                        if(AtomL==AtomK||AtomL==AtomI||AtomL<AtomK)continue;
                        add_imp({AtomI,AtomJ,AtomK,AtomL},Types,
                                Params,FoundCoords,Carts);    
                    }//Close AtomL imptorsion
                }//Close if imptorsion
            }//Close Atom K
            if(AtomJ<AtomI)continue; //Avoid 2x counting bond
            pair12.insert(std::make_pair(AtomI,AtomJ));
            add_bond({AtomI,AtomJ},Types,Params,FoundCoords);
        }//Close Atom J     
    }//Close AtomI
    for(size_t AtomI=0;AtomI<NAtoms;++AtomI){
        for(size_t AtomJ=AtomI+1;AtomJ<NAtoms;++AtomJ){
            auto pair=std::make_pair(AtomI,AtomJ);
            bool is12=pair12.count(pair),is13=pair13.count(pair),
                 is14=pair14.count(pair);
            if(is12||is13)continue;
            add_els({AtomI,AtomJ},Types,Params,FoundCoords,is14?chg14scale:1.0);
            add_lj({AtomI,AtomJ},Types,Params,FoundCoords,is14?vdw14scale:1.0);
        }
    }
    return FoundCoords;
}

void add_lj(const IVector& Atoms,const AtomTypes& Types,
             const ParamTypes& Params,CoordArray& FoundCoords,double scale){
    const IArray<2> ts={Types[Atoms[0]][0],Types[Atoms[1]][0]};
    const double sigmai=Params.at(LENNARD_JONES).at(sigma).at({ts[0]}),
          sigmaj=Params.at(LENNARD_JONES).at(sigma).at({ts[1]}),
          epsiloni=Params.at(LENNARD_JONES).at(epsilon).at({ts[0]}),
          epsilonj=Params.at(LENNARD_JONES).at(epsilon).at({ts[1]});
    const double avgsigma=sigmai+sigmaj,avgepsilon=sqrt(epsiloni*epsilonj);
    FoundCoords[LENNARD_JONES]->add_coord(Atoms,sigma,avgsigma,epsilon,avgepsilon*scale);
}

void add_els(const IVector& Atoms,const AtomTypes& Types,
             const ParamTypes& Params,CoordArray& FoundCoords,double scale){
    const IArray<2> ts={Types[Atoms[0]][1],Types[Atoms[1]][1]};
    const double qi=Params.at(ELECTROSTATICS).at(q).at({ts[0]}),
          qj=Params.at(ELECTROSTATICS).at(q).at({ts[1]});
    
    FoundCoords[ELECTROSTATICS]->add_coord(Atoms,Param_t::q,qi*qj*scale);
}

/* Begin indirection galore*/

void add_bond(const IVector& Atoms,const AtomTypes& Types,
              const ParamTypes& Params,CoordArray& FoundCoords){
    const IArray<2> ts={Types[Atoms[0]][0],Types[Atoms[1]][0]};
    const IVector ff_types={min(ts[0],ts[1]),max(ts[0],ts[1])};
    const double K_in=Params.at(BOND).at(K).at(ff_types);
    const double r0_in=Params.at(BOND).at(r0).at(ff_types);
    FoundCoords[BOND]->add_coord(Atoms,Param_t::K,K_in,Param_t::r0,r0_in);
}

void add_angle(const IVector& Atoms,
               const AtomTypes& Types,const ParamTypes& Params,
               CoordArray& FoundCoords){
    const IArray<3> ts={Types[Atoms[0]][0],Types[Atoms[1]][0],
                        Types[Atoms[2]][0]};
    const IVector ff_types={min(ts[0],ts[2]),ts[1],max(ts[0],ts[2])};
    const double K_in=Params.at(ANGLE).at(K).at(ff_types);
    const double r0_in=Params.at(ANGLE).at(r0).at(ff_types);
    FoundCoords[ANGLE]->add_coord(Atoms,K,K_in,r0,r0_in);
}

void add_torsion(const IVector& Atoms,
                 const AtomTypes& Types,const ParamTypes& Params,
                 CoordArray& FoundCoords){
    const IArray<4> ts={Types[Atoms[0]][0],Types[Atoms[1]][0],
                        Types[Atoms[2]][0],Types[Atoms[3]][0]};
    const bool already_ordered=ts[1]<ts[2]||((ts[1]==ts[2])&&ts[0]<ts[3]);
    const IVector ff_types=(already_ordered?IVector({ts[0],ts[1],ts[2],ts[3]}):
                                            IVector({ts[3],ts[2],ts[1],ts[0]}));
    double V_in=Params.at(TORSION).at(amp).at(ff_types);
    double n_in=Params.at(TORSION).at(n).at(ff_types);
    double phi_in=Params.at(TORSION).at(phi).at(ff_types);
    FoundCoords[TORSION]->add_coord(Atoms,amp,V_in,n,n_in,phi,phi_in);
    if(Params.at(TORSION).at(amp2).count(ff_types)){
        V_in=Params.at(TORSION).at(amp2).at(ff_types);
        n_in=Params.at(TORSION).at(n2).at(ff_types);
        phi_in=Params.at(TORSION).at(phi2).at(ff_types);
        FoundCoords[TORSION]->add_coord(Atoms,amp,V_in,n,n_in,phi,phi_in);
    }
    if(Params.at(TORSION).at(amp3).count(ff_types)){
        V_in=Params.at(TORSION).at(amp3).at(ff_types);
        n_in=Params.at(TORSION).at(n3).at(ff_types);
        phi_in=Params.at(TORSION).at(phi3).at(ff_types);
        FoundCoords[TORSION]->add_coord(Atoms,amp,V_in,n,n_in,phi,phi_in);
    }
}

//Code factorization for ordering the atoms in the improper torsions
template<typename T,typename T2,typename T3>
IArray<3> order_imp(size_t common,size_t v1,
           size_t v2,const T& spokes, const T2& mags, const T3& ts){
        const double angle1=std::acos(
            dot(spokes[common],spokes[v1])/(mags[common]*mags[v1]));
        const double angle2=std::acos(
            dot(spokes[common],spokes[v2])/(mags[common]*mags[v2]));
        const IArray<2> deg=angle1<angle2? IArray<2>({v1,v2}):IArray<2>({v2,v1});
        return ts[common]>ts[v1]? 
                 IArray<3>({deg[0],deg[1],common}):IArray<3>({common,deg[0],deg[1]});
}

using DArray=std::array<double,3>;

void add_imp(const IVector& Atoms,const AtomTypes& Types,
             const ParamTypes& Params,CoordArray& FoundCoords,const DVector& carts){
    const int tj=Types[Atoms[1]][0];//The type of the central atom
    //Eventually we need these sorted, but doing so destroys the map to the atoms
    IArray<3> ts={Types[Atoms[0]][0],Types[Atoms[2]][0],Types[Atoms[3]][0]};

    const std::array<DArray,3> spokes={//Distance from each atom to center
        diff(&carts[Atoms[1]*3],&carts[Atoms[0]*3]),
        diff(&carts[Atoms[1]*3],&carts[Atoms[2]*3]),
        diff(&carts[Atoms[1]*3],&carts[Atoms[3]*3])};
    const DArray mags={mag(spokes[0]),mag(spokes[1]),mag(spokes[2])};
    
    //Sort the atoms according to criteria described in manual
    IArray<3> final_order={0,1,2};//This will be the order of the periphereal atoms
    
    //Worry about cases where one atom type is degenerate
    if(ts[0]==ts[1])final_order=order_imp(2,0,1,spokes,mags,ts);
    else if(ts[0]==ts[2])final_order=order_imp(1,0,2,spokes,mags,ts);
    else if(ts[1]==ts[2])final_order=order_imp(0,1,2,spokes,mags,ts);
    else //Case where all unique (all same hits condition 1) 
        std::sort(final_order.begin(),final_order.end(),[&](size_t i, size_t j){
            return ts[i]<ts[j];}//End lambda
        );
    //Parameters are always looked up sorted
    std::sort(ts.begin(),ts.end());
    const IVector ff_types={ts[0],tj,ts[1],ts[2]};
    
    //Indices in final_order are only for periphereal atoms
    const IVector TempAtoms={Atoms[0],Atoms[2],Atoms[3]};
    const IVector FinalAtoms={TempAtoms[final_order[0]],Atoms[1],
                              TempAtoms[final_order[1]],TempAtoms[final_order[2]]};
    const double V_in=Params.at(IMPTORSION).at(amp).at(ff_types),
                 n_in=Params.at(IMPTORSION).at(n).at(ff_types),
               phi_in=Params.at(IMPTORSION).at(phi).at(ff_types);
    FoundCoords[IMPTORSION]->add_coord(FinalAtoms,amp,V_in,n,n_in,phi,phi_in);
}

} //End namespace FManII

