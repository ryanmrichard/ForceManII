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
namespace detail{class FFImpl;}

/** \brief  A struct to hold the details about a force field
 *
 * This class holds the definitions of a force field not a
 * force field applied to a molecule.  See [Force Fields in FManII](@ref ffdef)
 * for more details.
 *
 * \todo Fully hide implementation in pimpl member
 *
 */
class ForceField{
public:
    ///Type of a function that can order atoms
    typedef std::vector<size_t>(*orderer)(const std::vector<size_t>&);

    ///Type of a function that can combine parameters
    typedef double(*combiner)(const std::vector<double>&);

    using PTerm_t=std::pair<std::string,std::string>;

    ///Makes a force field that uses \p wildcard as the value for a wild_card
    ForceField(size_t wildcard=0);

    ///Frees memory associated with pimpl_
    ~ForceField();

    ///Deep copies other
    ForceField(const ForceField& other);

    ///Assigns this to a deep copy of other
    const ForceField& operator=(ForceField other);

    /** \brief Establishes that two force field terms are linked
     *
     *  Take for example the fact that we consider the van Der Waals 1-4
     *  interactions to be different than those of all other pairs.  Typically,
     *  one does not define parameters for say the 1-5 pairs, the 1-6 pairs, ...
     *  instead parameters are shared among these pairs.  We call this linking
     *  terms.  More specifically, if term 1 is linked to term 2, then when a
     *  user requests a parameter for term 1, the FF will first check parameters
     *  regestered under term 1, if it does not find it there it will then look
     *  for the parameter within term 2's set.
     *
     *  \param[in] term1 The parameter set to link
     *  \param[in] term2 The parameter set term1 is linked to
     */
    void link_terms(const FFTerm_t& term1,const FFTerm_t& term2);

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
     *  \param[in] atom_numbers A vector of the ordered sets of atom numbers in each coordinate
     *  \param[in] atom2type A mapping from atom number to atom type
     *  \param[in] skip_missing If true missing parameters will be ignored
     *
     *  \return The requested set of parameters
     */
    Vector assign_param(const FFTerm_t& term_type,
                        const std::string& parmi,
                        const std::vector<IVector>& atom_numbers,
                        const IVector& atom2type,
                        bool skip_missing)const;

    ///Checks for exact equality of all members
    bool operator==(const ForceField& other)const;

    ///Checks for inequality of any member
    bool operator!=(const ForceField& other)const{return !(*this==other);}

private:
    std::unique_ptr<detail::FFImpl> pimpl_;
};

///Functions for combining parameters
///@{
inline double mean(const Vector& params){
    return std::accumulate(params.begin(),params.end(),0.0)/params.size();
}

inline double product(const Vector& params){
    return std::accumulate(params.begin(),params.end(),1.0,std::multiplies<double>());
}

inline double geometric(const Vector& params){
    return std::pow(product(params),1.0/params.size());
}
///@}

///Some pre-defined order-er functions
///@{
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
///@}

} //End namespace FManII

