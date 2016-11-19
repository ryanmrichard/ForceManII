Terminology Used Throughout ForceManII                            {#terminology}
======================================

The purpose of this page is to clarify what we mean by some basic terms
throughout the manual and throughout the code.  I have tried to be consistent
but let me know if you find exceptions.

- atom number : Atoms in an input are listed in some order, whichever atom is
  listed first is said to have atom number 1, the second is atom number 2, etc.
- VDW type: I don't know how standard this is, and I'm willing to change the
  term if someone tells me the standard term for this.  Forcefields usually map
  an atom in a particular chemical enviornment, say an \f$sp^2\f$ hybridized 
  carbon in a carbonyl, to a single number, called the VDW type.  The VDW type
  is then mapped to another number called...
- atom type: This is less specific than the VDW type, say all \f$sp^2\f$ atoms
  map to the same atom type.  Properties that are more geometric in nature,
  bonds and angles, tend to list parameters by this type, whereas properties
  that are more chemical in nature, Coulomb and VDW potentials, tend to use the
  VDW type to list parameters.  If it helps clarify, atom types are only used
  internally and are read from the force field, the value you provide to
  ForceManII is always the VDW type

