Torsion Angle Class                                                   {#torsion}
===================

\warning Energies computed by ForcemanII for the torsional and improper
torsional energetic components may differ from those of other programs.  Read
on to see why this is.

# Torsion

Assume that atom 1 is bonded to atom 2; atom 2 is bonded to atom3, and atom 3 is
bonded to atom 4.  If we look down the 2-3 bond the 1-2 and 2-3 bonds will make
an angle.  This angle is the torsion angle.  Mathematically this can be viewed
as the angle between plane comprised of atoms 1,2, and 3 and the plane comprised
of the atoms 2,3,4.  Let \f$\vec{r_{21}}\f$ be a vector
parallel to the 1-2 bond pointing from atom 2 to atom 1, \f$\vec{r_{23}}\f$ be a
vector parallel to the 2-3 bond pointing from atom 2 to atom 3, and 
\f$\vec{r_{34}}\f$ be a vector parallel to the 3-4 bond pointing from atom 3 to
atom 4.  \f$\vec{n_1}=\vec{r_{21}}\times\vec{r_{23}}\f$ is a vector normal to 
the plane containing atoms 1, 2, and 3 and 
\f$\vec{n_2}=\vec{r_{34}}\times\vec{r_{23}}\f$
is a vector normal to the plane containing atoms 2, 3, and 4.  The dihedral
angle, \f$\phi\f$ is then given by:
\f[
\cos\theta=\frac{\vec{n_1}\cdot\vec{n_2}}{n_1 n_2}.
\f]
Alternatively:
\f[
\vec{n_1}\times\vec{n_2}=\frac{n_1 n_2\sin\theta \vec{n_3}}{n_3}
\f]
where \f$\vec{n_3}\f$ is a unit vector perpindicular to \f$\vec{n_1}\f$ and 
\f$\vec{n_2}\f$ and with direction consistent with the right-hand rule.  Since,
both of our cross products involved \f$\vec{r_{23}}\f$, we know that 
\f$\vec{n_3}=\frac{\vec{r_{23}}}{r_{23}}\f$.  Given
that \f$\vec{n_3}\f$ is a unit vector we may multiply both sides of 
the previous equation from the right by it to get:
\f[
\sin\theta=\frac{\left(\vec{n_1}\times\vec{n_2}\right)\cdot\vec{r_{23}}}
                {n_1 n_2 r_{23}}
\f]
Using the fact that \f$\tan\theta=\frac{\sin\theta}{\cos\theta}\f$ we also have
\f[
\tan\theta=\frac{\left(\vec{n_1}\times\vec{n_2}\right)\cdot\vec{r_{23}}}
{r_{23}}
\f]

As far as symmetries of the torsion angle are concerned, we have one degree of
freedom, we can read the sequence forward or backward.  Presently, we choose
to read the sequence in which ever manner makes atoms 2 and atoms 3 show up
in number order (order in input).  Parameters for torsions are expected to be
given in an analogous manner (atom type 2 should be less then atom type 3).

## Improper Torsion

Related to the torsion angle, the improper torsion angle measures the planarity
around an sp\f$^2\f$ hybridized atom.  Instead of being a non-branching path
i.e. 1 is connected to 2, which is connected to 3, which is connected to 4, the
atoms in an improper torsion are all bonded to atom 2.  It should be noted that
there are other conventions for which atom the central atom is (usually that it
is atom 3), but we choose it to be atom 2 as this is naturally where it falls
when looping over connectivity.  

It is also worth noting a few properties of the improper dihedral angle.  Assume
for a quadruple of atoms \$\lbrace a_1,a_2,a_3,a_4\rbrace\f$, we always compute 
our dihedral angle as the angle between the planes formed from the
first three atoms, \f$a_1, a_2, a_3\f$,  and the last three atoms,
\f$a_2, a_3, a_4\f$.  If we swap the atom originally labeled as 2 and the atom 
originally labeled as 3, which is the differing conventions for where the 
central atom goes, the angle changes sign.  Swapping atom 1 and atom 4 produces 
the same sign change.  However, swapping atoms 1 and 3 produces a different 
angle than the original, and swapping 3 and 4 produces yet another unique angle.
This is a problem as there are three unique angles that can be generated just by
interchanging the order in which the atoms are numbered.  After a long search, 
the best information I can find is on 
[this website](http://chempedia.info/info/165388/).  Basically they
suggest that the order of the orbtial atoms is that of their input and the
central atom is given third for AMBER and first for CHARMM.  Again, our code
requests that the central atom is given second.  

In order to ensure we always
compute the same angle, regardless of the order in which the atoms are specified
in the input file, we compute the torsions based on the atom type; specifically,
the atom that has the lowest atom type is the first atom, atom two is always the
central atom, atom 3 is the atom with the next lowest atom type, and atom 4 is
the atom with the highest atom type.  In the event of a tie, we number the atoms
with the same atom types in a clockwise fashion around the central atom.  This
is most easily done by computing the angle between the vector going from the
central atom to the unique atom and the vector going from the central atom to
the non-unique atom.  Lowest angle is closest to the unique angle in a
clockwise sense