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
#include <iostream>
#include <algorithm>
#include <sstream>

namespace FManII {

//Case-insensitive string compare
inline bool s_comp(const std::string& lhs, const std::string& rhs){
    if(lhs.size()!=rhs.size())return false;
    for(size_t i=0;i<lhs.size();++i)
        if (std::toupper(lhs[i]) != std::toupper(rhs[i]))return false;
    return true;
}

//Code factorization for setting the combination rule
inline void set_mean(std::pair<std::string,std::string> p,const std::string& token, ForceField& ff){
    if(s_comp(token,"GEOMETRIC"))ff.combrules[p]=FManII::geometric;
    else if(s_comp(token,"ARITHMETIC"))ff.combrules[p]=FManII::mean;
}

//Splits a line into strings based on whitespace
std::vector<std::string> tokenize(const std::string& line){
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    while(iss){
        std::string token;
        iss>>token;
        tokens.push_back(token);
    }
    return tokens;
}

template<typename Fterm,typename Fxn_t>
void generic_fill(ForceField&ff,
                     FFTerm_t ffterm,
                     const std::vector<std::string>& tokens,
                     size_t size,
                     Fxn_t fxn,
                     const std::string defaulttype,
                     const std::vector<std::pair<std::string,Vector>>& ps){
    ff.orderrules[ffterm]=fxn;
    if(!ff.paramtypes.count(ffterm))ff.paramtypes[ffterm]=defaulttype;
    if(!ff.terms.count(ffterm))ff.terms.emplace(ffterm,std::move(Fterm()));
    IVector vals;
    for(size_t i=1;i<=size;++i)vals.push_back(stoi(tokens[i]));
    vals=fxn(vals);
    for(size_t i=0;i<ps.size();++i)
        ff.params.add_param(ffterm,ps[i].first,vals,ps[i].second);
}

//Fills in the bond section of the FF
inline void parse_bond(const std::vector<std::string>& tokens,
                       ForceField& ff,
                       double k2au,
                       double ang2au){
    FFTerm_t ffterm=std::make_pair(Model_t::HARMONICOSCILLATOR,IntCoord_t::BOND);
    generic_fill<HarmonicBond>(ff,ffterm,tokens,2,pair_order,TypeTypes_t::CLASS,
        {{Param_t::K,{2.0*k2au*stod(tokens[3])}},
         {Param_t::r0,{ang2au*stod(tokens[4])}}});
}

//Fills in the angle section of the FF
inline void parse_angle(const std::vector<std::string>& tokens,
                        ForceField& ff,
                        double kcalmol2au,
                        double deg2rad){
    FFTerm_t ffterm=std::make_pair(Model_t::HARMONICOSCILLATOR,IntCoord_t::ANGLE);
    generic_fill<HarmonicAngle>(ff,ffterm,tokens,3,angle_order,TypeTypes_t::CLASS,
    {{Param_t::K,{2.0*kcalmol2au*stod(tokens[4])}},{Param_t::r0,{deg2rad*stod(tokens[5])}}});
}

//Fills in the torsion section of the FF
inline void parse_torsion(const std::vector<std::string>& tokens,
                          ForceField& ff,
                          double kcalmol2au,
                          double deg2rad){
    FFTerm_t ffterm=std::make_pair(Model_t::FOURIERSERIES,IntCoord_t::TORSION);
    Vector v,phi,n;
    for(size_t ti=0;ti<3;++ti){
        if(tokens.size()<=6+3*ti)break;//Null string on end of tokens
        v.push_back(stod(tokens[5+3*ti])*kcalmol2au);
        phi.push_back(stod(tokens[6+3*ti])*deg2rad);
        n.push_back(stod(tokens[7+3*ti]));
    }
    generic_fill<FourierTorsion>(ff,ffterm,tokens,4,torsion_order,
       TypeTypes_t::CLASS,{{Param_t::amp,v},{Param_t::phi,phi},{Param_t::n,n}});
}

//Fills in the improper torsion section of the FF
inline void parse_imp(const std::vector<std::string>& tokens,
                      ForceField& ff,
                      double kcalmol2au,
                      double deg2rad){
    const size_t i=stoi(tokens[1]),j=stoi(tokens[2]),
            k=stoi(tokens[3]),l=stoi(tokens[4]);
    //Tinker puts central atom third
    const std::vector<size_t> types=imp_order({i,k,j,l});
    FFTerm_t ffterm=std::make_pair(Model_t::FOURIERSERIES,IntCoord_t::IMPTORSION);
    ff.orderrules[ffterm]=imp_order;
    ff.paramtypes[ffterm]=TypeTypes_t::CLASS;
    if(!ff.terms.count(ffterm))
        ff.terms.emplace(ffterm,std::move(FourierImproperTorsion()));
    ff.params.add_param(ffterm,Param_t::amp,types,{stod(tokens[5])*kcalmol2au});
    ff.params.add_param(ffterm,Param_t::phi,types,{stod(tokens[6])*deg2rad});
    ff.params.add_param(ffterm,Param_t::n,types,{stod(tokens[7])});
}

//Fills in the 6-12 section of the force field
inline void parse_lj(const std::vector<std::string>& tokens,
                     ForceField& ff,
                     double kcalmol2au,
                     double ang2au,
                     bool isRMin,
                     bool isradius){
    const std::vector<size_t> types({stoul(tokens[1])});
    const double scale=ang2au*(isradius?2.0:1.0)*(isRMin?1.0:std::pow(2.0,1.0/6.0));
    for(const auto& ci:{IntCoord_t::PAIR,IntCoord_t::PAIR14}){
        FFTerm_t ffterm=std::make_pair(Model_t::LENNARD_JONES,ci);
        if(!ff.paramtypes.count(ffterm))ff.paramtypes[ffterm]=TypeTypes_t::CLASS;
        ff.params.add_param(ffterm,Param_t::sigma,types,{stod(tokens[2])*scale});
        ff.params.add_param(ffterm,Param_t::epsilon,types,{stod(tokens[3])*kcalmol2au});
        if(!ff.terms.count(ffterm)&&ci==IntCoord_t::PAIR14)
            ff.terms.emplace(ffterm,std::move(LJ14()));
        else if(!ff.terms.count(ffterm))
            ff.terms.emplace(ffterm,std::move(LJPair()));
    }

}

inline void parse_chg(const std::vector<std::string>& tokens,
                      ForceField& ff){
    const std::vector<size_t> types({stoul(tokens[1])});
    for(const auto& ci:{IntCoord_t::PAIR,IntCoord_t::PAIR14}){
        FFTerm_t ffterm=std::make_pair(Model_t::ELECTROSTATICS,ci);
        ff.paramtypes[ffterm]=TypeTypes_t::TYPE;
        ff.params.add_param(ffterm,Param_t::q,types,{stod(tokens[2])});
        if(!ff.terms.count(ffterm)&&ci==IntCoord_t::PAIR14)
            ff.terms.emplace(ffterm,std::move(Electrostatics14()));
        else if(!ff.terms.count(ffterm))
            ff.terms.emplace(ffterm,std::move(ElectrostaticsPair()));
    }
}

ForceField parse_file(std::istream&& file,
                      double kcalmol2au,
                      double ang2au,
                      double deg2rad){
    using std::stoi;using std::stod;
    const double k2au=kcalmol2au/(ang2au*ang2au);
    ForceField ff;
    std::string line;
    bool isRMin=true,isradius=true;
    while(std::getline(file,line)){
        if(line.size()==1)continue;
        auto tokens=tokenize(line);
        if(s_comp(tokens[0],"radiusrule"))
            set_mean({Model_t::LENNARD_JONES,Param_t::sigma},tokens[1],ff);
        else if(s_comp(tokens[0],"torsionunit"))
            ff.scale_factors[std::make_pair(Model_t::FOURIERSERIES,IntCoord_t::TORSION)]=stod(tokens[1]);
        else if(s_comp(tokens[0],"imptorunit"))
            ff.scale_factors[std::make_pair(Model_t::FOURIERSERIES,IntCoord_t::IMPTORSION)]=stod(tokens[1]);
        else if(s_comp(tokens[0],"vdwindex")){
            const auto& type=s_comp(tokens[1],"type")?TypeTypes_t::TYPE:TypeTypes_t::CLASS;
            ff.paramtypes[std::make_pair(Model_t::LENNARD_JONES,IntCoord_t::PAIR14)]=type;
            ff.paramtypes[std::make_pair(Model_t::LENNARD_JONES,IntCoord_t::PAIR)]=type;
        }
        else if(s_comp(tokens[0],"radiustype"))
            isRMin=s_comp(tokens[1],"R-MIN");
        else if(s_comp(tokens[0],"radiussize"))
            isradius=s_comp(tokens[1],"RADIUS");
        else if(s_comp(tokens[0],"epsilonrule"))
            set_mean({Model_t::LENNARD_JONES,Param_t::epsilon},tokens[1],ff);
        else if(s_comp(tokens[0],"vdw-14-scale"))
            ff.scale_factors[std::make_pair(Model_t::LENNARD_JONES,IntCoord_t::PAIR14)]=1.0/stod(tokens[1]);
        else if(s_comp(tokens[0],"chg-14-scale"))
            ff.scale_factors[std::make_pair(Model_t::ELECTROSTATICS,IntCoord_t::PAIR14)]=1.0/stod(tokens[1]);
        else if(s_comp(tokens[0],"atom"))
            ff.type2class[stoi(tokens[1])]=stoi(tokens[2]);
        else if(s_comp(tokens[0],"bond"))
            parse_bond(tokens,ff,k2au,ang2au);
        else if(s_comp(tokens[0],"angle"))
            parse_angle(tokens,ff,kcalmol2au,deg2rad);
        else if(s_comp(tokens[0],"imptors"))
            parse_imp(tokens,ff,kcalmol2au,deg2rad);
        else if(s_comp(tokens[0],"torsion"))
            parse_torsion(tokens,ff,kcalmol2au,deg2rad);
        else if(s_comp(tokens[0],"charge"))
            parse_chg(tokens,ff);
        else if(s_comp(tokens[0],"vdw"))
            parse_lj(tokens,ff,kcalmol2au,ang2au,isRMin,isradius);
    }
    for(const auto& pi:{IntCoord_t::PAIR,IntCoord_t::PAIR14}){
        auto ffterm=std::make_pair(Model_t::ELECTROSTATICS,pi);
        auto pterm=std::make_pair(Model_t::ELECTROSTATICS,Param_t::q);
        if(ff.terms.count(ffterm)){
            ff.combrules[pterm]=FManII::product;
        }
    }

    return ff;
}

}//End namespace
