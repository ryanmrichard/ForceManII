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
#pragma once

#include "ForceManII/FManIIDefs.hpp"
#include "ForceManII/FFTerm.hpp"
#include "ForceManII/ParameterSet.hpp"
#include <algorithm>
#include <set>
#include <unordered_map>

///Namespace for all code associated with ForceManII
namespace FManII {


/** \brief  A struct to hold the details about a force field
 *
 * This class holds the definitions of a force field not a
 * force field applied to a molecule.  See [Force Fields in FManII](@ref ffdef)
 * for more details.
 *
 */
struct ForceField{
    ///Type of a function that can order atoms
    typedef std::vector<size_t>(*orderer)(const std::vector<size_t>&);

    ///Type of a function that can combine parameters
    typedef double(*combiner)(const std::vector<double>&);

    using PTerm_t=std::pair<std::string,std::string>;

    ///Value of a parameter matching any type or class
    size_t wild_card=0;

    ParameterSet params;///<The complete set of parameters
    std::unordered_map<size_t,size_t> type2class;///<Map of atom type 2 atom class
    std::map<FFTerm_t,FFTerm> terms;///< The terms in the force field
    std::map<FFTerm_t,orderer> orderrules;///<How the parameters are ordered
    std::map<FFTerm_t,std::string> paramtypes;///<Does the parameter use class or type
    std::map<PTerm_t,combiner> combrules;///< How to combine parameters
    std::map<FFTerm_t,double> scale_factors;///<Scale terms by how much?

    /** \brief Given a set of internal coordinates assigns parameters
     *
     *  \param[in] term_type The model and intcoordinate of the term
     *  \param[in] parmi The type of parameter
     *  \param[in] coord The found internal coordinates
     *  \param[in] atom2type A mapping from atom number to atom type
     *  \param[in] skip_missing If true missing parameters will be ignored
     *
     *  \return The requested set of parameters
     */
    Vector assign_param(const FFTerm_t& term_type,
                        const std::string& parmi,
                        const InternalCoordinates& coord,
                        const IVector& atom2type,
                        bool skip_missing)const;

    ///Checks for exact equality of all members
    bool operator==(const ForceField& other)const{
        return (type2class==other.type2class &&
                orderrules==other.orderrules &&
                paramtypes==other.paramtypes &&
                combrules==other.combrules &&
                scale_factors==other.scale_factors &&
                terms==other.terms &&
                params==other.params);
    }

    ///Checks for inequality of any member
    bool operator!=(const ForceField& other)const{
        return !(*this==other);
    }
};

inline double mean(const Vector& params){
    return std::accumulate(params.begin(),params.end(),0.0)/params.size();
}

inline double product(const Vector& params){
    return std::accumulate(params.begin(),params.end(),1.0,std::multiplies<double>());
}

inline double geometric(const Vector& params){
    return std::pow(product(params),1.0/params.size());
}

///Some pre-defined order-er functions
inline std::vector<size_t> pair_order(const std::vector<size_t>& atoms){
    return {std::min(atoms[0],atoms[1]),std::max(atoms[0],atoms[1])};
}

inline std::vector<size_t> angle_order(const std::vector<size_t>& atoms){
    return {std::min(atoms[0],atoms[2]),atoms[1],std::max(atoms[0],atoms[2])};
}

inline std::vector<size_t> torsion_order(const std::vector<size_t>& atoms){
    if(atoms[1]<atoms[2]||(atoms[1]==atoms[2] && atoms[0]<atoms[3]))
        return atoms;
    return {atoms[3],atoms[2],atoms[1],atoms[0]};
}

inline std::vector<size_t> imp_order(const std::vector<size_t>& atoms){
    std::vector<size_t> types={atoms[0],atoms[2],atoms[3]};
    std::sort(types.begin(),types.end());
    return {types[0],atoms[1],types[1],types[2]};
}

} //End namespace FManII

