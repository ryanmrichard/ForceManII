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

#include "ForceManII/ForceField.hpp"
#include "ForceManII/Common.hpp"
#include<algorithm>
#include<sstream>
#include<iterator>

using namespace std;

namespace FManII{

namespace detail{

struct FFImpl{
    const size_t wild_card;

    FFImpl(size_t wildcard):
        wild_card(wildcard){}

    map<FFTerm_t,FFTerm_t> links;
};

}//end namespace detail

ForceField::ForceField(size_t wild_card):
    pimpl_(make_unique<detail::FFImpl>(wild_card)){}

ForceField::~ForceField(){}

ForceField::ForceField(const ForceField& other):
    params(other.params),type2class(other.type2class),terms(other.terms),
    orderrules(other.orderrules),paramtypes(other.paramtypes),
    combrules(other.combrules),scale_factors(other.scale_factors),
    pimpl_(make_unique<detail::FFImpl>(*other.pimpl_)){}

const ForceField& ForceField::operator=(ForceField other)
{
    std::swap(pimpl_,other.pimpl_);
    return *this;
}

bool ForceField::operator==(const ForceField& other)const
{
    return (type2class==other.type2class &&
            orderrules==other.orderrules &&
            paramtypes==other.paramtypes &&
            combrules==other.combrules &&
            scale_factors==other.scale_factors &&
            terms==other.terms &&
            params==other.params&&
            pimpl_->wild_card==other.pimpl_->wild_card);
}

/*void ForceField::set_type_2_class(size_t type,size_t class_in){
    pimpl_->type2class[type]=class_in;
}*/

void ForceField::link_terms(const FFTerm_t& term1,const FFTerm_t& term2){
    pimpl_->links[term1]=term2;
}

//Simple function to print an informative message about missing parameters
inline void no_param_error(const string& parami,
                           const IVector& types,
                           bool use_class){
    stringstream ss;
    ss<<"No parameter "+parami+" for "+(use_class?"class":"type")+" ";
    for(auto ti:types)ss<<ti<<" ";
    ss<<endl;
    throw runtime_error(ss.str());
}

//Code when we have to consider a combination rule
inline Vector handle_combrule(const FFTerm_t& term_type,
                            const string& parami,
                            const IVector& types,
                            const ParameterSet& params,
                            const map<FFTerm_t,FFTerm_t>& links,
                            bool skip_missing,
                            bool use_class)
{
    Vector in_params;
    for(auto ti:types){
        const IVector tiv({ti});
        Vector vs=params.get_param(term_type,parami,tiv);
        if(vs.size()==1)
            in_params.push_back(vs[0]);
        else{
            if(links.count(term_type))
                vs=params.get_param(links.at(term_type),parami,tiv);
            if(vs.size()==1)
                in_params.push_back(vs[0]);
            else if(!skip_missing)
                no_param_error(parami,{ti},use_class);
        }
    }
    return in_params;
}

template<typename fxn>
inline Vector handle_other(const FFTerm_t& term_type,
                           const string& parami,
                           const IVector& types,
                           const ParameterSet& params,
                           const map<FFTerm_t,FFTerm_t>& links,
                           size_t wildcard,
                           fxn orderer)
{
    const size_t n=types.size();
    vector<FFTerm_t> range({term_type});
    Vector vs;
    if(links.count(term_type))
        range.push_back(links.at(term_type));
    for(auto tt: range){
        for(size_t i=1;i<n;++i){
            vector<bool> is_zero(n);
            fill(is_zero.begin(),is_zero.begin()+i,true);
            do{
                IVector temp(types);
                for(size_t j=0;j<n;++j)
                    if(is_zero[j])temp[j]=wildcard;
                vs=params.get_param(tt,parami,orderer.at(tt)(temp));
                if(vs.size())return vs;
            }while(prev_permutation(is_zero.begin(),is_zero.end()));
        }
    }
}


Vector ForceField::assign_param(const FFTerm_t& term_type,
                                const string& parami,
                                const vector<IVector>& atom_numbers,
                                const IVector& atom2type,
                                bool skip_missing)const{
    const bool use_class=paramtypes.at(term_type)==TypeTypes_t::CLASS;
    Vector rv;
    for(auto typei:atom_numbers){
        IVector types;
        transform(typei.begin(),typei.end(),back_inserter(types),
            [&](size_t t){t=atom2type[t];return use_class?type2class.at(t):t;}
        );
        if(orderrules.count(term_type))
            types=orderrules.at(term_type)(types);
        auto prule=make_pair(term_type.first,parami);
        if(combrules.count(prule)){//Has combining rule?
            const Vector in_params=handle_combrule(term_type,parami,types,params,
                                          pimpl_->links,skip_missing,use_class);
            rv.push_back(combrules.at(prule)(in_params));
            continue;//Move to next coordinate
        }
        const bool is_vector=(term_type.first==Terms_t::FS_TORSION.first&&
                              term_type.second==Terms_t::FS_TORSION.second);
        //Try just grabbing the parameter
        Vector vs=params.get_param(term_type,parami,types);
        if(!vs.size())//Try wildcards and try links
            vs=handle_other(term_type,parami,types,params,
                            pimpl_->links,pimpl_->wild_card,orderrules);
        if(!vs.size() && !skip_missing)
            no_param_error(parami,types,use_class);
        if(is_vector)
            for(size_t i=0;i<3;i++)
                rv.push_back(vs.size()>i?vs[i]:0.0);
        else{
            DEBUG_CHECK(vs.size()<=1,"Wasn't expecting a vector");
            rv.push_back(vs.size()==1?vs[0]:0.0);
        }
     }
    return rv;
}

}//End namespace FManII
