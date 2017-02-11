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
#include <sstream>
#include <cmath>
#include <iostream>
#include <iomanip>

inline void test_header(const std::string& Msg){
    const std::string stars(80,'*');
    std::cout<<stars<<std::endl;
    size_t padding=(size_t)std::floor((80-4-Msg.size())/2.0);
    std::string lpad(padding,' '),rpad(padding+Msg.size()%2,' ');
    std::cout<<"* "<<lpad<<Msg<<rpad<<" *"<<std::endl;
    std::cout<<stars<<std::endl;
}

inline void test_footer(){
    const std::string stars(80,'*'),Msg("All Tests Passed!!!!!!!");
    std::cout<<stars<<std::endl;
    size_t padding=(size_t)std::floor((80-4-Msg.size())/2.0);
    std::string lpad(padding,' '),rpad(padding+Msg.size()%2,' ');
    std::cout<<"* "<<lpad<<Msg<<rpad<<" *"<<std::endl;
    std::cout<<stars<<std::endl;
}

#define TEST_THROW(fxn,msg)\
do{bool threw=false;try{fxn;}catch(const std::runtime_error& error){threw=true;}\
if(!threw)throw std::runtime_error(msg);}while(0)

/** Function for testing if two values are equievelent
 *  \param[in] actual The value you computed
 *  \param[in] theory The value you expected
 *  \param[in] Msg    A brief description of the check you are performing
 */
template<typename T>
inline void test_value(const T& actual,
                       const T& theory,
                       const std::string &Msg,
                       bool print=true){
    if(actual!=theory)
        throw std::runtime_error("Values do not match");
    if(print)std::cout<<Msg<<" : passed"<<std::endl;
}

/** Function for testing if two doubles are equivalent
 *  
 *  \param[in] actual The value you computed
 *  \param[in] theory The value you expected
 *  \param[in] tol    How close should actual and theory minimally be?
 *  \param[in] Msg    A brief description of the check you are performing
 */
inline void test_value(double actual,double theory,
                       double tol,const std::string& Msg ,bool print=true){
    if(std::fabs(actual-theory)>tol){
        std::stringstream m;
        m<<Msg<<": "<<std::setprecision(9)<<actual<<" does not match "<<
           theory<<" to within "<<tol<<std::endl;
        throw std::runtime_error(m.str());
    }
    if(print)std::cout<<Msg<<": passed"<<std::endl;
    
}

inline void compare_vectors(const std::vector<double>& actual,
                            const std::vector<double>& theory,
                            double tol,const std::string& msg){
    if(actual.size()!=theory.size())throw std::runtime_error("sizes do not match");
    for(size_t i=0;i<actual.size();++i)test_value(actual[i],theory[i],tol,msg,false);
    std::cout<<msg<<": passed"<<std::endl;    
}
