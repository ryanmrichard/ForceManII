Tinker Format                                                    {#tinkerformat}
=============

Finding the actual criteria for the format of the Tiker Parameter file is
somewhat difficult so I list it here.  The file itself is seperated into three
parts.  The first part is metadata about the force-field such as citations,
combining rules, scale factors, etc.  The second part lists all recognized atom
types and vdw types.  The third part lists the parameters arranged by atom or
vdw type.

\note As implemented in ForceManII, the Tinker Format is case-sensitive

## Force Field Metadata

This is the important metadata for a force field:

- `radiusrule` : can be either `ARITHMETIC` or `GEOMETRIC`.  This is the type of
   the average used to combine the van Der Waals radii
- `radiussize` : is the radius a `RADIUS` or a `DIAMETER`
- `epsilonrul` : same as `radiusrule` except for the well-depth
- `vdw-14-scale` : this is a scale factor that will be applied to the van Der 
   Waals interactions between 1-4 bonded atoms, i.e. ends of a torsion angle
- `chg-14-scale` : same as `vdw-14-scale` except for the charges
- `electric` : the electric perm                332.0522173
- `dielectric` : the dielectric constant (or absolute permitivity) for the force
   field in the same units as the electric permetivity

## Atom Types and VDW Types

This is the specification of the atom and vdw types the format is:

```
atom <vdw type> <atom type> <symbol> <description> <atomic number> <mass> <number of bonds
```

- `atom` : this is flag to specify that we are defining vdw and atom types
- `vdw type` : this is an integer specifying which vdw type we are defnining
- `atom type` : this is an integer specifying the atom type
- `symobl` : 1 to 2 characters for the atom type, is usually some mixture of the atomic
  symbol and a letter to describe the chemical environment like T for terminal
- `description` : a string giving a more readable description
- `atomic number` : The atomic number of the atom as it appears on the periodic
  table
- `mass` : the isotope averaged mass of the atom in Daltons, i.e. the mass you 
  find on the periodic table
- `number of bonds` : the usual number of bonds that an atom of this vdw type
  makes

## Parameters

The order of the parameter sections does not matter aside from the fact that if
a parameter occurs more than once, i.e. you specify the force constant for
atom types 1 and 2 twice, the last value will be used.

### Harmonic Bond Stretching Potential

```
bond   <atom type 1>  <atom type 2> <force constant> <r_0>
```

- `bond` : a flag to specify these parameters are for bonds
- `atom type 1` : the atom type of one of the two atoms
- `atom type 2` : the atom type for the other atom
- `force constant` : the force constant in kcal/(mol Angstrom^2), the 1/2 of
  the harmonic potential is included
- `r_0` : the equilibrium distance for the bond in Angstroms

### Angle-Bending Harmonic Potential Parameters

```
angle  <atom type 1> <vertex atom type> <atom type 2> <force constant> <theta_0>
```

- `angle` : a flag to specify that these are parameters for an angle
- `atom type 1` : one of the non-vertex atoms' types
- `vertex atom type` : the atom type of the vertex atom
- `atom type 2` : the other non-vertex atom's type
- `force constant` : the force constant in kcal/(mol*radians^2)
- `theta_0` : the equilibrium angle in degrees

### Torsion Fourier-Series Parameters
```
torsion <end 1> <center 1> <center 2> <end 2> <V> <gamma> <n> [<V> <gamma> <n> [<V> <gamma> <n>]] 
```

- `torsion` : a flag specifying these are parameters for a proper torsion
- `end 1` : the atom type of one of the ends of the torsion angle
- `center 1` : the atom type of the central atom bonded to the first end atom
- `center 2` : the atom type of the next central atom
- `end 2` : the atom type of the other end atom
- `V` : the amplitude of the Fourier Series component in kcal/mol
- `gamma` : the equlibrium value of the torsion in degrees
- `n` : the periodicity of the Fourier Series component
- Other values of `V`, `gamma`, and `n` can optionally be repeated up to another
  two times to specify additional terms in the Fourier Series for this 
  interaction

### Improper Torsion Fourier-Series Parameters
```
imptors <atom type 1> <atom type 2> <central atom type> <atom type 3> <V> <gamma> <n>
```

- `imptors` : a flag specifying these are parameters for an improper torsion
- `atom type 1` : the atom type of one of the orbital atoms
- `atom type 2` : the atom type of one of the other orbital atoms
- `central atom type` : the atom type of the central atom
- `atom type 3` : the atom type of the last orbital atom
- `V` : the amplitude of the Fourier Series component in kcal/mol
- `gamma` : the equlibrium value of the torsion in degrees
- `n` : the periodicity of the Fourier Series component

### Charge-Charge Potential Parameters

```
charge <vdw type> <q>
```

- `charge` : a flag specifying these parameters are for a charge
- `vdw type` : the vdw type of the atom
- `q` : the charge in atomic units, i.e. number of electrons


### Lennard-Jones 6-12 Potential Parameters

```
vdw           <atom_type>               <radius> <epsilon>
```

- `vdw` : a flag to specify that these are parameters for van Der Waals 
  interaction
- `atom_type` : the atom type this parameter is for
- `radius` : the radius (or diameter depending on the metadata provided) in
  Angstroms
- `epsilon` : the well depth in kcal/mol




