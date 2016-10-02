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
#include <ForceManII/HarmonicOscillator.hpp>
#include "TestMacros.hpp"
#include "ubiquitin.hpp"
#include <iostream>
#include <cmath>

int main(int argc, char** argv){
    std::cout<<"Testing Harmonic Oscillator force-field term"<<std::endl;
    FManII::HarmonicOscillator HO;
    
#ifndef NDEBUG
//Test debug checks
    try{
        std::vector<double> a(3,0.0),b(2,0.0),c;
        c=HO.deriv(a,b);
        throw std::runtime_error("Check that bond.size()==ks.size() failed");
    }
    catch (const std::runtime_error& error){}
#endif
    
    //Energy check
    std::vector<double> Energy=HO.deriv(ubiquitin_bonds,ubiquitin_ks);
    test_value(Energy[0],ubiquitinbond_e,1e-5,"Bond Energy");
    
    
    return 0;
} //End main

