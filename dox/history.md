History and Release Information                                       {#history}
===============================

ForceManII was originally a code written to replace Q-Chem's original MM code,
which is entitled ForceMan (short for force-field manager).  In 2012, I (Ryan M.
Richard) started work on a revised updated version which I termed ForceManII.
The predominent goal of that rewrite was to add Hessians so that the MM code
could be used as a guess for geometry optimizations.  This was work that I
started as a graduate student at Ohio State (where I was
part of the Q-Chem development team) and finished at the beginning of my
postdoc at Georgia Institue of Technology (where I became a Psi4 developer).  In
an effort to still maintain ForceManII, the decision was made to make it a stand
alone library that could be called from either Q-Chem or Psi4.  This is the
result of that rewrite.

From the begining the design considerations of ForceManII are:

- Designed to be driven not the driver
  - Most MM packages assume they are in command as do most QM packages.  This
    leads to a power struggle.  ForceManII is designed to be callable either
    from a more fully featured MM package or a QM package.
- Designed to be extensible
  - QM/MM computations are becoming exceedingly popular.  Newer QM/MM procedures
    usually impose small modifications to existing procedures, modifications
    that we aim to allow the user to enforce on us without modifying the code
  - New force fields or reparameterizations of old ones are common we aim to be
    able to accomodate thes new developments without having to modify the core
    of ForceManII
  - A hot topic in QM methods development is on-the-fly force field
    parameterization and more physically motivated force fields we strive to
    provide a scaffolding to faciliate this.
- Designed to use modern object-oriented paradigms
  - Automation of much of the messy logic plauging older MM codes
  - Inheritance to reduce code duplication
  - STL usage for faster/better algorithms


## Versioning

# Release -1.0 (circa 2013)

The initial ForceManII project is propietary code in Q-Chem.  It includes some
nice features such as analytic Hessians, a more object-oriented design,
(attempted) means of easily extending, and reading parameters from files.  Aside
from the ideas of this release nothing remains of it.

# Release 1.0 (currently in pre-alpha)

This is the first release of ForceManII actually intended for use.  The main
features are:

- Support for several force fields
  - AMBER99
  - CHARMM22
  - OPLS-AA
- C++11 API that allows a high-level of flexibility on the user's part
  - Can use custom force fields and parameters
    - Either by manually filling in a `ForceField` class or by providing such a
      force field in an std::istream formatted according to Tinker's .prm
      conventions
    - Relatedly one may manually tweak the assigned parameters in the `ParamSet`
      instance before they are actually used
  - Ability to manually define internal coordinates
    - User may provide their own filled `CoordArray` instance
  - Disjoint implementations of model potentials and internal coordinates
    - Adding a new bond-stretching potential (*e.g.* a fifth order Taylor
      series) involvies only adding a new potential.  The derivatives coded for
      the distance between the two points (the internal coordinate) are reusable
      and need not be recoded.  This also works in reverse; to add a new
      internal coordinate, one need not rederive the potential used for it.
- Minimal input (only the Cartesian coordinates of the molecule, the atom types
  of each atom, and the connectivity are required)
  - Future releases aim to reduce this to just the Cartesian coordinates (and
    the atomic numbers)
  - Avoids the need for things like topology files, which are unfamiliar to QM
    people and replaces it with quantities that are commonly found in most QM
    packages (plus atom types)
- Analytic gradients

