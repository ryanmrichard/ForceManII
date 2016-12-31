Force Field Definition in ForceManII        {#ffdef}
====================================

ForceManII breaks from the traditional ways of thinking about many MM related
ideas.  What is a force field is one of them.  This page is meant to explain my
logic and also introduce you to the ForceField class.

## Abstract Background

Overall, in MM we think of our energy, \f$E\f$, as being a function of a set of
Cartesian coordinates, the vector of which we denote \f$\vec{x}\f$.  \f$E\f$ is
then thought of as being composed of some small integer number of terms,
\f$n_{\text{terms}}\f$, such that
\f$E(\vec{x})=\sum_{=1}^{n_\text{terms}}E_i(\vec{x})\f$, where \f$E_i\f$ is the
\f$i\f$-th term.  In general, the terms do not directly depend on \f$\vec{x}\f$,
but rather on the mapping of \f$\vec{x}\f$ to some other vector \f$\vec{q}\f$.
If \f$\vec{x}\f$ is the vector of Cartesian coordinates then \f$\vec{q}\f$ is
a vector of the internal coordinates (bonds, angles, etc.).  It is somewhat
helpful to partition \f$\vec{q}\f$ into segments where each segment is all
internal coordinates of one type (say bonds).  Thus the mapping to the
\f$i\f$-th type of internal coordinate is \f$\vec{q}_i(\vec{x})\f$.

Of course the hallmark of a force field is its set of emperical parameters,
which we label \f$\vec{t}\f$.  In general there will be one set of parameters
for each interal coordinate type, the one associated with the \f$i\f$-th
internal type being \f$\vec{t}_i\f$.  Although somewhat simple to do in
practice, abstractly defining the mapping between the system and
the parameters is no small feat, in part due to a lack of a standard in the MM
community.  In the most general case, one first maps \f$\vec{x}\f$ to a set of
chemical moieties (functional groups, amino acids, DNA bases, etc.).  Each of
the chemical moieties is then mapped to an integer.  We have no need of the
intermediate mapping presently and define \f$n(\vec{x}\f$ to be the mapping of
the system to a set of integers, with each integer being affiliated with a
particular chemical moiety.  The integers are the atom types.  This mapping is
usually too much to parameterize, so similar chemical moieties are parameterized
in groups.  This leads to a second map \f$N[n(\vec{x})]\f$, which generates what
we call our chemcial classes, but in general \f$N\f$ is to be thought of as the
mapping to the final integer.  These classes are for atoms, not internal
coordinates so
we need another mapping that defines associates these classes with internal
coordinates; conceptually this is similar to \f$\vec{q}\f$ so we define the
mapping to the \f$i\f$-th internal coordinate as \f$\vec{q}_i(N[n(\vec{x})])\f$.
Now this is mapped to the final set of parameters, which gives rise to the
mapping: \f$\vec{t}_i[\vec{q}_i(N[n(\vec{x})])]\f$.  Unfortunately, this isn't
quite complete as sometimes, depending on the internal coordinate type, the
parameters in \f$\vec{t}_i\f$ are then mapped to another set of parameters,
\f$\vec{t}^\prime_i\f$.


It is thus more appropriate to write our energy expansion as
\f[
E(\vec{x};\mathbf{T},\vec{c})=
\sum_{=1}^{n_\text{terms}}c_iE_i(
\mathbf{q}(\vec{x});\mathbf{t}^\prime(\mathbf{t}[\mathbf{q}(N[n(\vec{x})])])).
\f]
where we have also allowed for a reweighting of the terms so that energy term
\f$i\f$, has a weight \f$c_i\f$ and introduced bold face quantities to be the
set of that quantity for all internal coordinates.  The definition of all of the
above quantities and mappings defines a force field.  Like I said, I think this
is a rather unique perspective on what a force field is, but it is one that
falls out of the computer science.

## Parts of a Force Field

From the discussion above we conclude that the input to a force field is a
molecular system and the output is an energy or energy derivative.  The most
general force field includes definitions of the required mappings, as well
as the parameters and weights.  This means we need you to define:

- A set of energy terms
  - The form of these terms are usually simple: harmonic oscillator,
    Fourier Series, etc.
  - The same form may be used for multiple terms, but with a different
    set of parameters (*e.g* harmonic bond-stretching vs. harmonic
    angle-bending)
- A set of weights for each term
  - Usually these are unity, unless the term is for a pair of atoms seperated
    by four or five bonds.
- An initial set of parameters for each internal coordinate
  - Usually a map from the union of the atom type/classes to some scalars, but
    may also be to vectors
- A means of mapping the initial set of parameters to the final set
  - Can be an identity mapping
  - Usually only used for combining epsilons and sigmas of Lennard Jones
    potentials so that one need not store all pair-wise combinations of them
  - Usually is something simiple like arthemetic or geometric mean
  - Note this requires a mapping of the atomic charges to their product in the
    case of the Coulomb term.  This is because to work correctly in this
    framework all energy terms take the parameters for the internal coordinate
    and not the individual atoms.
- A list of which terms use atom type and which use atom class
- A mapping of atom type to atom class

ForceManII will (can) take care of:
- Computing the internal coordinates (given the system and connectivity)
- Mapping the internal coordinates to the final integer (given the ff, inteneral
  coordinates, and types of each atom)
- Mapping the final integers to the final parameters (given the ff)
- Computing the derivative (given the internal coordinates and final parameters)

In short, other than defining a force field (we provide some of the more
popular ones for you) you need to give us the Cartesian coordinates,
connectivity, and atom types of your system to use the library.
