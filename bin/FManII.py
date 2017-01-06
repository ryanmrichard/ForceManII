"""This script includes things common to the other included Python scripts
"""

import math
import os
import sys
#Constants as defined in ForceManII
ang2au=1.889725989
j2cal=4.184
kcalmol2au=1/627.5096
deg2rad=math.pi/180
k2au=kcalmol2au/(ang2au**2)
anglek2au=kcalmol2au

#Recognized param types
params={"k":"FManII::Param_t::K",
        "r0":"FManII::Param_t::r0",
        "v":"FManII::Param_t::amp",
        "phi":"FManII::Param_t::phi",
        "n":"FManII::Param_t::n",
        "q":"FManII::Param_t::q",
        "sigma":"FManII::Param_t::sigma",
        "epsilon":"FManII::Param_t::epsilon"
}

#Recognized internal coordinate types
intcoords={"bond":"FManII::IntCoord_t::BOND",
       "pair13":"FManII::IntCoord_t::PAIR13",
       "pair14":"FManII::IntCoord_t::PAIR14",
       "pair":"FManII::IntCoord_t::PAIR",
       "angle":"FManII::IntCoord_t::ANGLE",
       "torsion":"FManII::IntCoord_t::TORSION",
       "imp":"FManII::IntCoord_t::IMPTORSION",
}

#Recognized model types
models={"ho":"FManII::Model_t::HARMONICOSCILLATOR",
        "fs":"FManII::Model_t::FOURIERSERIES",
        "cl":"FManII::Model_t::ELECTROSTATICS",
        "lj":"FManII::Model_t::LENNARD_JONES",
}

#Recognized ways of combining parameters
comb_rules={"ARITHMETIC":"FManII::mean",
           "GEOMETRIC":"FManII::geometric",
           "PRODUCT":"FManII::product"}


#Recognized atom type mappings
typetypes={"type":"FManII::TypeTypes_t::TYPE",
           "class":"FManII::TypeTypes_t::CLASS"
           }

#Recognized force field terms
ffterms={"hb":(models["ho"],intcoords["bond"]),
        "ha":(models["ho"],intcoords["angle"]),
        "hi":(models["ho"],intcoords["imp"]),
        "ft":(models["fs"],intcoords["torsion"]),
        "fi":(models["fs"],intcoords["imp"]),
        "cl14":(models["cl"],intcoords["pair14"]),
        "cl":(models["cl"],intcoords["pair"]),
        "lj14":(models["lj"],intcoords["pair14"]),
        "lj":(models["lj"],intcoords["pair"]),
        "ub":(models["ho"],intcoords["pair13"])}

def check(cond,msg):
    if not cond: raise RuntimeError(msg)

def cross(v1,v2):
    """Defines cross product so we don't depend on numpy"""
    return [v1[1]*v2[2]-v1[2]*v2[1],
            v1[2]*v2[0]-v1[0]*v2[2],
            v1[0]*v2[1]-v1[1]*v2[0]]

def dot(v1,v2):
    """Defines dot product so we don't depend on numpy"""
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]

def diff(v1,v2):
    """Defines the difference between two vectors"""
    return [v1[k]-v2[k] for k in range(0,3)]

def mag(v1):
    """Magnitude of a vector"""
    return math.sqrt(sum(v1[k]**2 for k in range(0,3))) 


def read_sys(xyz_file):
    """Given the path to a Tinker-style .xyz file reads the coordinates,
       connectivity, and parameter types into lists.  These lists are then
       returned"""
    carts=[];connect=[];param_num=[]
    with open(xyz_file,"r") as f:
        next(f)
        for line in f:
            tokenized=line.split()
            carts.append([float(i)*ang2au for i in tokenized[2:5]])
            param_num.append(int(tokenized[5]))
            connect.append([int(i)-1 for i in tokenized[6:]])#Tinker starts at 1
    return carts,connect,param_num
