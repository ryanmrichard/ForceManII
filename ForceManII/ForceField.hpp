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
/** \file ForceField.hpp
 * 
 * \version 0.1
 * \date November 27, 2016 at 11:59 AM (EST)
 *  
 * Original Author: \author Ryan M. Richard (ryanmrichard1<at>gmail.com)
 * 
 * Additional contributions by:
 *
 */

#ifndef FMANII_FORCEFIELD_HPP
#define FMANII_FORCEFIELD_HPP

#include "ForceManII/FManIIDefs.hpp"
#include <set>

///Namespace for all code associated with ForceManII
namespace FManII {

/** \brief  A struct to hold the details about a force field
 *
 * This class holds the definitions of a force field not a
 * force filed applied to a molecule.
 */
struct ForceField{
    ParamTypes params;///<The complete set of parameters
    std::set<IntCoord_t> terms;///< The types of terms
    double chg14scale=1.0;///<What to scale 1-4 chg-chg by
    double vdw14scale=1.0;///<What to scale 1-4 vdw by
    std::map<Param_t,CombRule_t> combrules;
};

} //End namespace FManII

#endif /* End header guard */

