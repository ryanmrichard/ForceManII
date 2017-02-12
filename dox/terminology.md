Terminology Used Throughout ForceManII                            {#terminology}
======================================

The purpose of this page is to clarify what we mean by some basic terms
throughout the manual and throughout the code.  I have tried to be consistent
but let me know if you find exceptions.

- atom number : Atoms in an input are listed in some order, whichever atom is
  listed first is said to have atom number 0, the second is atom number 1, etc.
- atom type: I don't know how standard this and the next terms are, but this
  is the language used by Tinker.  Forcefields usually map
  an atom in a particular chemical enviornment, say an \f$sp^2\f$ hybridized 
  carbon in a carbonyl, to a single number, called the atom type.  The atom type
  is then mapped to another number called...
- class type: This is less specific than the atom type, say all \f$sp^2\f$ atoms
  map to the same class type.  Properties that are more geometric in nature,
  bonds and angles, tend to list parameters by this type, whereas the atom type
  appears to only be used for Coulomb potentials.  If it helps clarify, class
  types are only used internally and are read from the force field, the value
  you provide to ForceManII is always the atom type
- force field term: This is a contribution to the energy caused by, for example,
  stretching a bond, bending an angle, etc.
- Cartesian coordinate: The XYZ coordinates of an atom and the standard way of
  specifying a molecule in electronic structure theory
- Connectivity table: A table of the atoms an atom is bonded to
- internal coordinate: This is the length of a bond, the angle of an angle, the
  angle of a torsion angle, etc. The internal coordinates depend on the
  Cartesian coordinates of each atom and the bond assignments
- model potential: This is the functional form of the energy term.  It is
  usually something like a harmonic oscillator or Fourier series.

