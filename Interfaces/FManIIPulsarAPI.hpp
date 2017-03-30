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

#include<pulsar/modulebase/EnergyMethod.hpp>
#include<pulsar/modulebase/PropertyCalculator.hpp>

namespace FManII {


class FFPulsar:public pulsar::EnergyMethod{
public:
    FFPulsar(ID_t id): EnergyMethod(id){}
    pulsar::DerivReturnType deriv_(size_t Order,const pulsar::Wavefunction& wfn);
};

class FFTermPulsar:public pulsar::EnergyMethod{
public:
    FFTermPulsar(ID_t id): EnergyMethod(id){}
    pulsar::DerivReturnType deriv_(size_t Order,const pulsar::Wavefunction& wfn);
};

class FFCharges:public pulsar::PropertyCalculator{
public:
    FFCharges(ID_t id): PropertyCalculator(id){}
    std::vector<double> calculate_(unsigned int deriv,
                                   const pulsar::Wavefunction & wfn);
};

}//End namespace FManII
