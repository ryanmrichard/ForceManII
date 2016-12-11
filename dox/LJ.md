Lennard-Jones (6-12) Potential                                             {#LJ}
==============================

## Theory

The van Der Waals term of a force-field accounts for the exchange repulsion of
the atoms as well as the dispersive attraction of the atoms.  The repulsive
term is modeled as being proportional to \f$\vec{r}^{-12}\f$, where 
\f$\vec{r}\f$ is the distance between the atoms and the attractive part is
modeled after dispersion which is known to behave as \f$\vec{r}^{-6}\f$.  This
means the final form is:
\f[
E=\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-
                \left(\frac{\sigma}{r}\right)^{6}\right]
\f]
All forms of this potential contain two parameters, \f$\epsilon\f$, which is
the well-depth, and \f$\sigma\f$ the distance between the two atoms for which
the interaction is zero.  This form was introduced by John Lennard-Jones, hence
the name.  For whatever reason there is usually a 4 out front giving:
\f[
E=4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12}-
                \left(\frac{\sigma}{r}\right)^{6}\right]
\f]
Sometimes \f$\sigma\f$ is actually given as the distance between the two atoms
for which the interaction is minimized.  Taking the first derivative of the
previous equation and setting it to 0:
\f[
0=\frac{-24\epsilon}{r}\left[\left(\frac{\sigma}{r}\right)^{12}-
      \left(\frac{\sigma}{r}\right)^{6}\right]
\f]
fromw which we deduce that the minimum distance is:
\f[
r=2^{1/6}\sigma
\f]
hence \f$\sigma\f$ at the minimum, \f$\sigma^\prime\f$, is \f$2^{-1/6}\sigma\f$.
Rewriting the potential in terms of \f$sigma^\prime\f$:
\f[
E=\epsilon\left[\left(\frac{\sigma^\prime}{r}\right)^{12}-
                2*\left(\frac{\sigma^\prime}{r}\right)^{6}\right]
\f]


Force fields usually list \f$\epsilon\f$ and \f$\sigma\f$ individually for each
atom and then use a combining rule to to form the final value.  The combining
rule is either the arthmetic or geometric mean.  One tricky point is whether or
not the value listed for \f$\sigma\f$ is given as the radius or the diameter of
the van Der Waals sphere around it; it needs to be a diameter.  To see this,
consider \f$\sigma^\prime\f$ for a homodiatomic molecule, the ideal distance is
then twice the radius of \f$\sigma^\prime\f$.  For the herterodiatomic we are
assuming the ideal distance is an average of the two homodiatomic distances,
*i.e.* the averages of the diameters.

## Conventions

The following are our conventions for a Lennard-Jones term:
- \f$\sigma^\prime\f$ is used instead of \f$\sigma\f$
- \f$\sigma^\prime\f$ is given as a diameter
- \f$\sigma^\prime\f$ and \f$\epsilon\f$ are pre-averaged before being passed
  to the Lennard-Jones class
- As usual all input is expected to be in atomic units
