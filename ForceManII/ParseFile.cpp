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
inline void set_mean(Param_t p,const std::string& token, ForceField& ff){
    if(s_comp(token,"GEOMETRIC"))ff.combrules[p]=GEOMETRIC;
    else if(s_comp(token,"ARITHMETIC"))ff.combrules[p]=ARITHMETIC;
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

//Fills in the bond section of the FF
inline void parse_bond(const std::vector<std::string>& tokens,
                       ParamTypes& params,
                       double k2au,
                       double ang2au){
    const size_t i_=stoi(tokens[1]),j_=stoi(tokens[2]);
    const size_t i=std::min(i_,j_),j=std::max(i_,j_);
    const std::vector<size_t> types({i,j});
    params[BOND][K][types]=2.0*stod(tokens[3])*k2au;
    params[BOND][r0][types]=stod(tokens[4])*ang2au;
}

//Fills in the angle section of the FF
inline void parse_angle(const std::vector<std::string>& tokens,
                        ParamTypes& params,
                        double kcalmol2au,
                        double deg2rad){
    const size_t i_=stoi(tokens[1]),j=stoi(tokens[2]),k_=stoi(tokens[3]);
    const size_t i=std::min(i_,k_),k=std::max(i_,k_);
    const std::vector<size_t> types({i,j,k});
    params[ANGLE][K][types]=stod(tokens[4])*kcalmol2au;
    params[ANGLE][r0][types]=stod(tokens[5])*deg2rad;
}

//Fills in the improper torsion section of the FF
inline void parse_imp(const std::vector<std::string>& tokens,
                      ParamTypes& params,
                      double kcalmol2au,
                      double deg2rad){
    const size_t i_=stoi(tokens[1]),j_=stoi(tokens[2]),
            k=stoi(tokens[3]),l_=stoi(tokens[4]);
    std::vector<size_t> types_={i_,j_,l_};
    std::sort(types_.begin(),types_.end());
    const std::vector<size_t> types({types_[0],types_[1],k,types_[2]});
    params[IMPTORSION][amp][types]=stod(tokens[5])*kcalmol2au;
    params[IMPTORSION][phi][types]=stod(tokens[6])*deg2rad;
    params[IMPTORSION][n][types]=stod(tokens[7]);
}

//Fills in the torsion section of the FF
inline void parse_torsion(const std::vector<std::string>& tokens,
                          ParamTypes& params,
                          double kcalmol2au,
                          double deg2rad){
    size_t i=stoi(tokens[1]),j=stoi(tokens[2]),k=stoi(tokens[3]), l=stoi(tokens[4]);
    if(j<k){
        k=j;
        j=stoi(tokens[3]);
        i=l;
        l=stoi(tokens[1]);
    }
    const std::vector<size_t> types({i,j,k,l});
    params[TORSION][amp][types]=stod(tokens[5])*kcalmol2au;
    params[TORSION][phi][types]=stod(tokens[6])*deg2rad;
    params[TORSION][n][types]=stod(tokens[7]);
    if(tokens.size()>8){
        params[TORSION][amp2][types]=stod(tokens[8])*kcalmol2au;
        params[TORSION][phi2][types]=stod(tokens[9])*deg2rad;
        params[TORSION][n2][types]=stod(tokens[10]);
        if(tokens.size()>11){
            params[TORSION][amp3][types]=stod(tokens[8])*kcalmol2au;
            params[TORSION][phi3][types]=stod(tokens[9])*deg2rad;
            params[TORSION][n3][types]=stod(tokens[10]);
        }
    }
}

//Fills in the 6-12 section of the force field
inline void parse_lj(const std::vector<std::string>& tokens,
                     ParamTypes& params,
                     double kcalmol2au,
                     double ang2au,
                     bool isRMin,
                     bool isradius){
    const std::vector<size_t> types({stoi(tokens[1])});
    const double scale=ang2au*(isradius?2.0:1.0)*(isRMin?1.0:std::pow(2.0,1.0/6.0));
    params[LENNARD_JONES][sigma][types]=stod(tokens[2])*scale;
    params[LENNARD_JONES][epsilon][types]=stod(tokens[3])*kcalmol2au;
}

ForceField parse_file(std::istream&& file,
                      double kcalmol2au,
                      double ang2au,
                      double deg2rad){
    using std::stoi;using std::stod;
    const double k2au=kcalmol2au/(ang2au*ang2au);
    ForceField ff;
    ParamTypes params;
    std::map<size_t,size_t> type2class;
    std::string line;
    bool isRMin=true,isradius=true;
    while(std::getline(file,line)){
        if(line.size()==0)continue;
        auto tokens=tokenize(line);
        if(s_comp(tokens[0],"raidusrule"))
            set_mean(sigma,tokens[1],ff);
        else if(s_comp(tokens[0],"radiustype"))
            isRMin=s_comp(tokens[1],"R-MIN");
        else if(s_comp(tokens[0],"radiussize"))
            isradius=s_comp(tokens[1],"RADIUS");
        else if(s_comp(tokens[0],"epsilonrule"))
            set_mean(epsilon,tokens[1],ff);
        else if(s_comp(tokens[0],"vdw-14-scale"))
            ff.vdw14scale=1.0/stod(tokens[1]);
        else if(s_comp(tokens[0],"chg=14-scale"))
            ff.chg14scale=1.0/stod(tokens[1]);
        else if(s_comp(tokens[0],"atom"))
            type2class[stoi(tokens[1])]=stoi(tokens[2]);
        else if(s_comp(tokens[0],"bond"))
            parse_bond(tokens,params,k2au,ang2au);
        else if(s_comp(tokens[0],"angle"))
            parse_angle(tokens,params,kcalmol2au,deg2rad);
        else if(s_comp(tokens[0],"imptors"))
            parse_imp(tokens,params,kcalmol2au,deg2rad);
        else if(s_comp(tokens[0],"torsion"))
            parse_torsion(tokens,params,kcalmol2au,deg2rad);
        else if(s_comp(tokens[0],"charge"))
            params[ELECTROSTATICS][q][{stoi(tokens[1])}]=stod(tokens[2]);
        else if(s_comp(tokens[0],"vdw"))
            parse_lj(tokens,params,kcalmol2au,ang2au,isRMin,isradius);
    }
    ff.params=params;
    for(auto type: params)
      ff.terms.insert(type.first);
    return ff;
}

}//End namespace
