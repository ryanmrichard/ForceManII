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

namespace FManII{

namespace detail{
    ///Class that will implement the details of the ParameterSet class
    class ParameterSetImpl;
    ///Class that will implement the details of the ParameterSetItr class
    class ParameterSetItrImpl;
}

/** \brief This is the API to the class that stores our parameters
 */
class ParameterSet{

    ///Class for looping over all parameters as if they were tuples
    struct ParameterSetItr{

        ParameterSetItr();///<Does nothing

        ///Makes an iterator that points to start of host if begin==true
        ParameterSetItr(const ParameterSet* host,bool begin);

        ///Deep copies other
        ParameterSetItr(const ParameterSetItr& other);

        ///Assigns this to a deep copy of other
        const ParameterSetItr& operator=(ParameterSetItr other);

        ///Frees memory
        ~ParameterSetItr();

        ///Type of the tuple this returns
        using value_type=std::tuple<FFTerm_t,std::string,IVector,Vector>;

        ///Get the tuple currently pointed at
        value_type operator*()const;

        ///Increments and then returns iterator
        ParameterSetItr& operator++(int);

        ///Increments after returning iterator
        ParameterSetItr operator++();

        ///True if this iterator is the same as other
        bool operator==(const ParameterSetItr& other)const;

        ///True if the iterators are not the same
        bool operator!=(const ParameterSetItr& other)const{
            return !(*this==other);
        }

    private:
        ///The actual implementation of the ParameterSet class
        std::unique_ptr<detail::ParameterSetItrImpl> pimpl_;
    };
public:
    ///Type of a constant iterator to this class
    using const_iterator=ParameterSetItr;

    const_iterator begin()const;///Returns an iterator to the first parameter
    const_iterator end()const;///Returns an iterator past the last parameter

    ParameterSet();///<Initializes the pimpl instance
    ~ParameterSet();///<Releases unique_ptr
    ParameterSet(const ParameterSet& other);///<Deep copies the pimpl instance

    ///copy/swap assignment operator
    const ParameterSet& operator=(ParameterSet other);

    /** \brief Adds a parameter to the set
     *
     *  \param[in] term The type of the ff term these parameters are for
     *                  e.g. harmonic-bond
     *  \param[in] type What sort of parameter is this (e.g. force constant)
     *  \param[in] atoms What atom types or classes is this affiliated with?
     *  \param[in] params The parameters
     *
     */
    void add_param(const FFTerm_t& term,
                   const std::string& type,
                   const IVector& atoms,
                   const Vector& params);

    ///Overload to reduce compile time for scalar parameters
    void add_param(const FFTerm_t& term,
                   const std::string& type,
                   const IVector& atoms,
                   double params);
    ///Overload to reduce compile time for scalar parameters and types
    void add_param(const FFTerm_t& term,
                   const std::string& type,
                   size_t atoms,
                   double params);

    /** \brief Gets a parameter accounting for the fact that the requested
     *  parameter may involve wildcards
     *
     *  \param[in] term The type of the ff term these parameters are for
     *                  e.g. harmonic-bond
     *  \param[in] type What sort of parameter is this (e.g. force constant)
     *  \param[in] atoms What atom types or classes is this affiliated with?
     *  \return The parameters as a vector or an empty vector if the parameter
     *  is not in the set
     *
     */
    Vector get_param(const FFTerm_t& term,
                     const std::string& type,
                     const IVector& atoms)const;

    ///Returns true if two sets of parameters are equal
    bool operator==(const ParameterSet& other)const;

    ///Returns true if two sets of parameters are not equal
    bool operator!=(const ParameterSet& other)const{
        return !(*this==other);
    }

private:
    std::unique_ptr<detail::ParameterSetImpl> pimpl_;
};




}//End namespace
