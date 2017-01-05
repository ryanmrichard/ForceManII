/*
 * Copyright (C) 2016 Ryan M. Richard.
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
#include "ForceManII/ModelPotential.hpp"

///Namespace for all code associated with ForceManII
namespace FManII {

class FFTerm{
    std::shared_ptr<const ModelPotential> model_;///< The model to use
public:
    ///Makes a new FF term given the model and coordinates it depends on
    FFTerm(std::shared_ptr<const ModelPotential> m,
           const std::vector<std::string>& cs):
        model_(std::move(m)),coords(cs){}

    ///Returns the model
    const ModelPotential& model()const{return *model_;}

    ///Returns the key for the term
    FFTerm_t name()const{return {model_->name,coords[0]};}

    const std::vector<std::string> coords;///<Coordinates this term depends on

    ///Given the parameters for this term and a set of internal coordinates
    ///computes the derivative
    Vector deriv(size_t order,const std::map<std::string,Vector>& ps,
                              const CoordArray& cs)const{
        std::vector<Vector> incoords;
        for(auto ci:coords)incoords.push_back(cs.at(ci)->get_coords());
        return model_->deriv(order,ps,incoords);
    }

    ///True if both terms have the same model and type of coordinates
    bool operator==(const FFTerm& other)const{
        return(*model_==*other.model_ &&
               coords==other.coords);
    }

    ///True if either term has a different model or type of coordinates
    bool operator!=(const FFTerm& other)const{
        return !this->operator==(other);
    }
};

} //End namespace FManII


