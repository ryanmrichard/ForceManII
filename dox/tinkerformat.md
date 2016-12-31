Tinker Format                                                    {#tinkerformat}
=============

Tinker has invented two of its own file formats: .prm and .xyz (the tinker .xyz
is not quite the same as the more prevelant standard version seen in
electronic structure theory).  This page explains the two file formats and
represents what we expect them to look like at the moment (there's nothing
stopping Tinker from modifying them after all).

- [Tinker Parameter File](#file_tinker_param)
- [Tinker .xyz File](#file_tinker_xyz)


## The Tinker Parameter File {#file_tinker_param}

Finding the actual criteria for the format of the Tiker parameter file is
somewhat difficult.  What is listed here is reverse engineered from the
parameter files included with Tinker.  The file itself is seperated into three
parts.  The first part is metadata about the force-field such as citations,
combining rules, scale factors, etc.  The second part lists all recognized atom
types and class types.  The third part lists the parameters arranged by atom or
class type. A line starting with `#` character is considered a comment line.
`!!` is also considered a comment character, but is typically used at the end of
the line.

The ForceManII API provides the function `parse_file` which given an
`std::istream` instance (the base class of either an input file or an input
stringstream) will parse the Tinker formated stream and provide you a
ForceField instance.  The parser follows the rules listed in this document.

\note As implemented in ForceManII, the Tinker Format is case-insensitive

### Force Field Metadata

This is the important metadata for a force field:

- `radiusrule` : can be either `ARITHMETIC` or `GEOMETRIC`.  This is the type of
   the average used to combine the van Der Waals radii
- `radiussize` : is sigma a `RADIUS` or a `DIAMETER`
- `radiustype` : is sigma for the minimum `R-MIN` or the zero?
- `epsilonrul` : same as `radiusrule` except for the well-depth
- `vdw-14-scale` : this is a scale factor that will be applied to the van Der 
   Waals interactions between 1-4 bonded atoms, i.e. ends of a torsion angle.
   The actual value applied will be 1.0 over this value
- `chg-14-scale` : same as `vdw-14-scale` except for the charges
- `electric` : the electric perm                332.0522173
- `dielectric` : the dielectric constant (or absolute permitivity) for the force
   field in the same units as the electric permetivity

### Atom Types and Class Types

This is the specification of the atom and class types the format is:

```
atom <atom type> <class type> <symbol> <description> <atomic number> <mass>
<number of bonds>
```

- `atom` : this is flag to specify that we are defining vdw and atom types
- `atom type` : this is an integer specifying which atom type we are defnining
- `class type` : this is an integer specifying the class type
- `symbol` : 1 to 2 characters describing the atom type, is usually some mixture
   of the atomic
  symbol and a letter to describe the chemical environment like T for terminal
- `description` : a string giving a more readable description
- `atomic number` : The atomic number of the atom as it appears on the periodic
  table
- `mass` : the isotope averaged mass of the atom in Daltons, i.e. the mass you 
  find on the periodic table
- `number of bonds` : the usual number of bonds that an atom of this atom type
  makes

### Parameters

The order of the parameter sections does not matter aside from the fact that if
a parameter occurs more than once, i.e. you specify the force constant for
class types 1 and 2 twice, the last value will be used.

#### Harmonic Bond Stretching Potential

```
bond   <class type 1>  <class type 2> <force constant> <r_0>
```

- `bond` : a flag to specify these parameters are for bonds
- `class type 1` : the class type of one of the two atoms
- `class type 2` : the class type for the other atom
- `force constant` : the force constant in kcal/(mol Angstrom\f$^2\f$), the 1/2
  of the harmonic potential is included
- `r_0` : the equilibrium distance for the bond in Angstroms

#### Angle-Bending Harmonic Potential Parameters

```
angle  <class type 1> <vertex class type> <class type 2> <force constant> <theta_0>
```

- `angle` : a flag to specify that these are parameters for an angle
- `class type 1` : one of the non-vertex atoms' class types
- `vertex atom type` : the atom type of the vertex atom
- `class type 2` : the other non-vertex atom's class type
- `force constant` : the force constant in kcal/(mol*radians^2)
- `theta_0` : the equilibrium angle in degrees

#### Torsion Fourier-Series Parameters
```
torsion <end 1> <center 1> <center 2> <end 2> <V> <gamma> <n> [<V> <gamma> <n> [<V> <gamma> <n>]] 
```

- `torsion` : a flag specifying these are parameters for a proper torsion
- `end 1` : the class type of one of the ends of the torsion angle
- `center 1` : the class type of the central atom bonded to the first end atom
- `center 2` : the class type of the next central atom
- `end 2` : the class type of the other end atom
- `V` : the amplitude of the Fourier Series component in kcal/mol
- `gamma` : the equlibrium value of the torsion in degrees
- `n` : the periodicity of the Fourier Series component
- Other values of `V`, `gamma`, and `n` can optionally be repeated up to another
  two times to specify additional terms in the Fourier Series for this 
  interaction

#### Improper Torsion Fourier-Series Parameters
```
imptors <class type 1> <class type 2> <central class type> <class type 3> <V> <gamma> <n>
```

- `imptors` : a flag specifying these are parameters for an improper torsion
- `atom type 1` : the atom type of one of the orbital atoms
- `atom type 2` : the atom type of one of the other orbital atoms
- `central atom type` : the atom type of the central atom
- `atom type 3` : the atom type of the last orbital atom
- `V` : the amplitude of the Fourier Series component in kcal/mol
- `gamma` : the equlibrium value of the torsion in degrees
- `n` : the periodicity of the Fourier Series component

#### Charge-Charge Potential Parameters

```
charge <vdw type> <q>
```

- `charge` : a flag specifying these parameters are for a charge
- `vdw type` : the vdw type of the atom
- `q` : the charge in atomic units, i.e. number of electrons


#### Lennard-Jones 6-12 Potential Parameters

```
vdw           <atom_type>               <radius> <epsilon>
```

- `vdw` : a flag to specify that these are parameters for van Der Waals 
  interaction
- `atom_type` : the atom type this parameter is for
- `radius` : the radius (or diameter depending on the metadata provided) in
  Angstroms
- `epsilon` : the well depth in kcal/mol

## The Tinker .xyz File {#file_tinker_xyz}


