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
#include "ForceManII/ParameterSet.hpp"
#include "ForceManII/Common.hpp"
#include <unordered_map>
#include <limits>

template <class T>
inline void hash_combine(size_t& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
namespace std
{
    template<typename T>
    struct hash<vector<T>>
    {
        typedef vector<T> argument_type;
        typedef std::size_t result_type;
        result_type operator()(argument_type const& in) const
        {
            size_t size = in.size();
            size_t seed = 0;
            for (size_t i = 0; i < size; i++)
                //Combine the hash of the current vector with the hashes of the previous ones
                hash_combine(seed, in[i]);
            return seed;
        }
    };
}

using namespace std;

namespace FManII{

//Map between an ordered tuple to its parameter value (value may be a tensor)
using IndexedParam=unordered_map<IVector,Vector>;

//Map between a type of parameter and the indexed parameters
using Type2Index_t=unordered_map<string,IndexedParam>;

//Full map of a ff's parameters
using Term2Type_t=map<FFTerm_t,Type2Index_t>;

using value_type=tuple<FFTerm_t,string,IVector,Vector>;

namespace detail{

struct ParameterSetImpl{
    Term2Type_t params;
};

struct ParameterSetItrImpl{
    using OuterItr=typename Term2Type_t::const_iterator;
    OuterItr outer_itr;
    using MiddleItr=typename Type2Index_t::const_iterator;
    MiddleItr middle_itr;
    using InnerItr=typename IndexedParam::const_iterator;
    InnerItr inner_itr;
    const ParameterSetImpl* ps;

    void reset(){
        auto itr1=ps->params.begin();
        bool is_good=itr1!=ps->params.end();
        middle_itr=(is_good?itr1->second.begin():MiddleItr());
        is_good=is_good && itr1->second.begin()!=itr1->second.end();
        inner_itr=(is_good?itr1->second.begin()->second.begin():InnerItr());
    }

    ParameterSetItrImpl(const ParameterSetImpl* psin,bool begin):
        outer_itr(begin? psin->params.begin():psin->params.end()),
        ps(psin){
        reset();
    }

    void next(){
        if(++inner_itr!=middle_itr->second.end())return;
        else if(++middle_itr!=outer_itr->second.end()){
            inner_itr=middle_itr->second.begin();
            return;
        }
        else{
            ++outer_itr;
            if(outer_itr!=ps->params.end()){
                middle_itr=outer_itr->second.begin();
                inner_itr=middle_itr->second.begin();
            }
            else reset();
            return;
        }
    }

    value_type operator*()const{
        return make_tuple(outer_itr->first,
                          middle_itr->first,
                          inner_itr->first,
                          inner_itr->second);
    }

    bool operator==(const ParameterSetItrImpl& other)const{
        return (ps==other.ps &&
                outer_itr==other.outer_itr &&
                middle_itr==other.middle_itr &&
                inner_itr==other.inner_itr
                );
    }
};

}//End namespace detail

using Itr_t=typename ParameterSet::const_iterator;

Itr_t::ParameterSetItr(const ParameterSet* host,bool begin):
    pimpl_(make_unique<detail::ParameterSetItrImpl>(&(*host->pimpl_),begin)){}

Itr_t::ParameterSetItr(){}

Itr_t::~ParameterSetItr(){}

Itr_t& Itr_t::operator++(int){
    pimpl_->next();
    return *this;
}

Itr_t Itr_t::operator++(){
    auto rv=Itr_t(*this);
    pimpl_->next();
    return rv;
}

value_type Itr_t::operator*()const{
    return pimpl_->operator*();
}

Itr_t::ParameterSetItr(const Itr_t& other):
    pimpl_(make_unique<detail::ParameterSetItrImpl>(*other.pimpl_)){}

bool Itr_t::operator==(const Itr_t& other)const{
    return *pimpl_==*other.pimpl_;
}

const Itr_t& Itr_t::operator=(Itr_t other){
    std::swap(this->pimpl_,other.pimpl_);
    return *this;
}

Itr_t ParameterSet::begin()const{return Itr_t(this,true);}
Itr_t ParameterSet::end()const{return Itr_t(this,false);}

ParameterSet::ParameterSet():pimpl_(make_unique<detail::ParameterSetImpl>()){}

ParameterSet::~ParameterSet(){}

ParameterSet::ParameterSet(const ParameterSet& other):
    pimpl_(make_unique<detail::ParameterSetImpl>(*other.pimpl_)){}

const ParameterSet& ParameterSet::operator=(ParameterSet other){
    std::swap(this->pimpl_,other.pimpl_);
    return *this;
}
bool ParameterSet::operator==(const ParameterSet& other)const{
    return pimpl_->params==other.pimpl_->params;
}

void ParameterSet::add_param(const FFTerm_t& term,const string& type,
                             const IVector& atoms,double params)
{
    add_param(term,type,atoms,Vector(1,params));
}

void ParameterSet::add_param(const FFTerm_t& term,const string& type,
                             size_t atoms,double params)
{
    add_param(term,type,IVector(1,atoms),Vector(1,params));
}

void ParameterSet::add_param(const FFTerm_t& term,const string& type,
                             const IVector& atoms,const Vector& params)
{
    pimpl_->params[term][type][atoms]=params;
}

Vector ParameterSet::get_param(const FFTerm_t& term, const string& type,
                               const IVector& atoms)const{
    if(pimpl_->params.count(term))
        if(pimpl_->params.at(term).count(type))
            if(pimpl_->params.at(term).at(type).count(atoms))
                return pimpl_->params.at(term).at(type).at(atoms);
    return {};
}

}//End namespace FManII
