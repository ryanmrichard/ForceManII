Force Field Terms in FManII                                 {#force_field_terms}
===========================

The heart and soul of ForceManII is the FFTerm class.  This class is what turns
a model potential and a set of internal coordinates into an energy or an
energy derivative.  Extending ForceManII thus comes down to being able to write
new model potentials and possibly new internal coordinates.  This page aims to
elucidate the details of the class and hopefully provide you enough information
to design your own FFTerms if you feel so inclined.

## Theoretical Background

Let's start with some basic theory.  For our present purposes, each term in a
force field can be thought of as a function that returns an energy given an
internal coordinate, \f$Q\f$ (for the moment, we ignore cross terms that depend
on multiple internal coordinates).  We usually think of the coordinate as a
single number, the distance, the angle, etc.  However, when we want derivatives
of this energy we want them with respect to Cartesian coordinates.
Consequentially it helps to think of this internal coordinate as a function,
that when given a vector of Cartesian coordinates, \f$\vec{q}\f$, returns a
number.  Regardless, calling our energy function \f$E\f$ the first derivative of
the energy with respect to \f$\vec{q}\f$ is readily computed with the chain
rule:
\f[
    \frac{\partial E[Q(\vec{q})]}{\partial \vec{q}}=
         \frac{\partial E[Q(\vec{q})]}{\partial Q(\vec{q})}
         \frac{\partial Q(\vec{q})}{\partial \vec{q}},
\f]

The second derivative with respect to \f$\vec{q}\f$ follows from the product
rule:
\f[
    \frac{\partial^2 E[Q(\vec{q})]}{\partial\vec{q}^2}=
    \frac{\partial^2 E[Q(\vec{q})]}{\partial Q(\vec{q})^2}
    \frac{\partial Q(\vec{q})}{\partial \vec{q}}+
    \frac{\partial E[Q(\vec{q})]}{\partial Q(\vec{q})}
    \frac{\partial Q(\vec{q})^2}{\partial \vec{q}^2}
\f]
and we ignore higher order derivatives.

What this mathematical model tells us is the final derivate really depends on
two things, the derivative of the energy's functional form and the derivative of
the coordinate.  In fact, we should be able to code these up separtely. This
also facilitates code reuse because some energy functions are used
in multiple force field terms (*e.g.* it is not uncommon for both the bond and
the angle terms to be harmonic oscillators) and some internal coordinates are
associated with multiple force field terms (*e.g.* both the bond stretching and
the Lennard-Jones terms depend on the distance between two points).

## Class Design

\TODO Describe why members go where they are and why they are the types they are


## Extending ForceManII

In order to extend ForceManII's capabilities new FFTerms have to be defined and
then compiled into your source code.  For example, pretending ForceManII does
not know how to compute bond-stretching via a harmonic oscillator:

~~~{.cpp}
#include <ForceManII/FManII.hpp> //Main API of ForceManII

//Define a model potential that implements a harmonic oscillator by inheritance
class HarmonicOscillator : public FManII::ModelPotential
{
public:
    //Define a default constructor that passes a list of parameter names to the
    //base classs.  Parameter names will be namespace protected by the potential
    //so generic names are fine.  In this example "K" will be our force constant
    //and "R0" will be the equilibrium bond length
    HarmonicOscillator():
        ModelPotential({Param_t::k,Param_t::r0}){}

    //Implment the deriv function (Vector is a typedef of std::vector<double>)
    FManII::Vector deriv(size_t order,
                         const std::map<FManII::Param_t,FManII::Vector>& in_params,
                         const std::vector<FManII::Vector>& in_coords)
    {
        //Normally this would be implemented in a .cpp file, but for brevity we
        //are doing it in the header and ignoring all energy derivatives

        //One element long vector initialized to 0 that will be our energy
        std::vector<double> return_value(1,0.0);
        for(size_t i=0; i<in_coords[0].size();++i)
        {
            const double qi=(in_coords[0][i]-in_params.at(Param_t::r0)[i]);
            return_value[0]+=0.5*in_params.at(Param_t::k)[i]*qi*qi;
        }

        return return_value;
    }//End deriv
};//End Harmonic Oscillator class
~~~

\TODO Describe how to derive from the InternalCoordinates class
