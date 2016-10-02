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
/** \file TestMacros.hpp
 * 
 * \version 0.1
 * \date October 2, 2016 at 2:20 PM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FMANII_TESTMACROS_HPP
#define FMANII_TESTMACROS_HPP
#include <sstream>
#include <cmath>

/** Function for testing if two doubles are equivalent
 *  
 *  \param[in] actual The value you computed
 *  \param[in] theory The value you expected
 *  \param[in] tol    How close should actual and theory minimally be?
 *  \param[in] Msg    A brief description of the check you are performing
 */
inline void test_value(double actual,double theory,
                       double tol,const std::string& Msg ){
    if(std::fabs(actual-theory)>tol){
        std::stringstream m;
        m<<Msg<<": "<<actual<<" does not match "<<theory<<" to within "
         <<tol<<std::endl;
        throw std::runtime_error(m.str());
    }
}

#endif /* End header guard */

