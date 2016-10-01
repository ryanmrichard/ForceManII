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
/** \file Common.hpp
 * 
 * \version 0.1
 * \date October 1, 2016 at 1:23 PM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FMANII_COMMON_HPP
#define FMANII_COMMON_HPP
#include <stdexcept>

///Namespace for all code associated with ForceManII
namespace FManII {

//Macro that checks if a condition is true, if it's not macro throws 
//(only for Debug builds)
#ifdef NDEBUG
#define DEBUG_CHECK(cond,msg)
#else
#define DEBUG_CHECK(cond,msg)do{if(!cond)throw std::runtime_error(msg);}while(0)
#endif /* End ndebug*/

} //End namespace FManII

#endif /* End header guard */

