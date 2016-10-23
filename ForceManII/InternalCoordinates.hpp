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

/* 
 * File:   InternalCoordinates.hpp
 * Author: Ryan M. Richard <ryanmrichard1 at gmail.com>
 *
 * Created on October 9, 2016, 11:46 AM
 */

#ifndef INTERNALCOORDINATES_HPP
#define INTERNALCOORDINATES_HPP
#include "ForceManII/FManII.hpp"

namespace FManII {

/** \brief A base class designed to provide basic functionality to things like
 *         bonds, angles, etc.
 * 
 *  Energetic terms in a force field rely on determining how far an internal
 *  coordinate, or IntCoord for short, is from an equilibrium value.  These
 *  IntCoords are things like bonds, angles, torsions, etc. and the force field
 *  maps the energetic cost of displacing these IntCoords to simple functional
 *  forms.  Hence it is key that we have a class to describe these internal
 *  coordinates.
 * 
 *  Each instance of this class is designed to work with a single system.  If 
 *  you need to have IntCoords for multiple systems make multiple instances.
 *  Basic workflow with this class is give it a natoms by 3 array of Cartesian
 *  coordinates at creation.  At creation the order of the atoms is irrelevant,
 *  but after creation all further input with this class will refer to that
 *  ordering.  Next an algorithm above this class locates the internal
 *  coordinate of choice (i.e. using the connectivity you find say a bond) you
 *  then add that bond to this class like:
 *  \code{.cc}
 *  IntCoords Bonds(System);
 *  Bonds.add_coord({1,2},Param_t::K,value_of_K);
 *  \endcode
 *  and similarly for other coordinates.
 * 
 *  
 *  
 */
class IntCoords {   
protected:
    ///The number of atoms this coordinate depends on
    size_t natoms_;
    
    ///Typedef for the set of atoms in an IntCoord
    using Atoms_t=std::vector<size_t>;
    
    ///Typedef for a vector of doubles
    using VDouble=std::vector<double>;
    
    ///The cartesian coordinates of the system in a.u.
    std::shared_ptr<const VDouble> carts_;
    
    ///The actual length of the bond, or angle of the angle, etc.
    VDouble value_;
    
    ///The atom numbers associated with each value
    std::vector<Atoms_t> atoms_;
    
    ///The parameters associated with this set of internal coordinates
    std::map<Param_t,VDouble> params_;
    
    /** \brief Overridden by the derived class so that it computes the i-th 
     *         derivative of the coordinate which involves the specified atoms
     * 
     *   We provide you the derivative order and the coordinate number.  At the
     *   moment I am assuming that while computing the value some of the
     *   parameters can be baked into the value to simplify things later.
     *   Specifically I am assuming that the 0-th order derivative can just
     *   return the difference between the actual value and the equilibrium
     *   distance for bonds, angles, and Urey-Bradley terms.
     * 
     *   \param[in] deriv_i What order-th derivative should you return?
     *   \param[in] coord_i Which coordinate are you taking the derivative of?
     *   \return A \f$(3N)^{i}\f$ long array of the derivatives in the order:
     *   x of the first atom in the coord, y of the first atom in the coord, z
     *   of the first atom in the coord, x of the second atom in the coord,...
     *   along each of the \f$i\f$ dimensions of the derivative.
     */ 
    virtual VDouble compute_value(size_t deriv_i, Atoms_t coord_i)const=0;
    
    ///Recursion end point
    void add_coord(const Atoms_t& coord){
        value_.push_back(compute_value(0,coord)[0]);
    }
    
        ///\p Carts is coordinates of the molecule in a.u. in a NAtoms X 3 array
    IntCoords(size_t NAtoms,std::shared_ptr<const VDouble> Carts):
        natoms_(NAtoms),carts_(Carts){}
    
public:  
    
    /** \brief The call to add a coordinate to this IntCoord
     * 
     *  Don't be offput by the signature using this function is pretty easy.
     *  The first argument is the atoms in this coordinate. The remaining 2N
     *  arguments are the N parameters associated with this coordinate.  They
     *  must be in the order parameter type, value of that parameter, next
     *  parameter type, value of the next parameter, etc.  You will get a
     *  compile time error if you do not supply pairs for each parameter or if
     *  you mix up the order.
     * 
     *  Point is usage looks like this:
     *  \code{.cc}
        //Assume these are already made instances of IntCoords
        IntCoords MyBonds,MyAngles;
      
        //Adds a bond between atoms 1 and 2 with a force constant of 3.4
        MyBonds.add_coord({1,2},Param_t::K,3.4);
      
        //Adds an angle among atoms 1, 2, and 3 with force constant 3.3
        MyAngles.add_coord({1,2,3},Param_t::K,3.3);
        \endcode
     * 
     *  \param[in] coord An array of the atoms in the coordinate.  Numbers must
     *                   correspond to the order given at instance creation
     *  \param[in] param The type of the parameter, e.g. force constant
     *  \param[in] pvalue The value of the parameter in a.u.
     *  \param[in] other_params Remaining parameter type, parameter value pairs
     */
    template<typename...args>
    void add_coord(const Atoms_t& coord,Param_t param,
                   double pvalue,args...other_params){
        static_assert(!(sizeof...(args)%2),
                "Must provide pairs of parameters and values");
        params_[param].push_back(pvalue);
        add_coord(coord,other_params...);
    }
    
    ///Returns the array of coordinates
    const VDouble& values()const{return value_;}
    
    ///Returns the array of a particular parameter
    const VDouble& params(Param_t pt)const{return params_.at(pt);}
    
    ///Virtual Dtor so compiler doesn't complain, but does nothing
    virtual ~IntCoords(){}
    

};


}//End namespace FManII

#endif /* INTERNALCOORDINATES_HPP */

