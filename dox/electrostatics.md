Electrostatics                                                 {#electrostatics}
==============

\warning The included electrostatics term in ForceManII works a little different
than in most other programs. Read on to find out why.

\note At the moment I have only coded up charge-charge interactions according to
Coulomb's law.  I have left the name of this term more general anticipating the
desire to code up more electrostatic terms.

In classical electrostatics, given a field of \f$N\f$ point charges, such that
the \f$j\f$-th one is located at \f$\vec{r}_j\f$ and has charge \f$q_j\f$, the
energy to insert an additional point charge of charge \f$q_i\f$ at a position 
\f$\vec{r}_i\f$ is given by:
\f[
E_i=\sum_{j=0}^N \frac{q_iq_j}{|\vec{r}_i-\vec{r}_j|},
\f]
which is Coloumb's law for point charges (up to some constants that depend on
the unit choice).  The quantity:
\f[
V(\vec{r}_{i};q_j)=\sum_{j=0}^N \frac{q_j}{|\vec{r}_i-\vec{r}_j|}
\f]
is the electrostatic potential at the point \f$\vec{r}_i\f$.  Coulomb's law
gives us the energy to place the \f$i\f$-th charge in the current field.  If
we want the total energy to assemble the field it is given by considering the
placement of each charge and is then given by:
\f[
E=\sum_i E_i=\sum_{i<j}\frac{q_iq_j}{|\vec{r}_i-\vec{r}_j|}
\f]

In ForceManII all models are given a vector of the internal coordinates they
depend on (in this case the distance between two points) and the parameters
for each of these internal coordinates.  However, the form of Coulomb's law
gives the parameters (the charges) in terms of each atom, not a given distance.
We can get around this by defining a new parameter, \f$t_{ij}=q_iq_j\f$, for
each distance.  Consequentially, charge-charge interactions used as input in
ForceManII should be \f$t_{ij}\f$.  This can be handled in an automated fashion
by defining the combining rule of your force field's charge parameters to be
their product.  Thus our final energy expression is:
\f[
E(\vec{r};\vec{q})=\sum_{i<j}\frac{t_{ij}}{|\vec{r}_i-\vec{r}_j|}
\f]
the derivative with respect to \f$r_{ij}\f$ is:
\f[
\frac{\partial E}{\partial r_{ij}}=\frac{-t_{ij}}{|\vec{r}_i-\vec{r}_j|^2}
\f]
and the second derivative is:
\f[
\frac{\partial^2 E}{\partial r_{ij}^2}=\frac{2t_{ij}}{|\vec{r}_i-\vec{r}_j|^3}
\f]
the \f$n\f$-th derivative is:
\f[
\frac{\partial^n E}{\partial r_{ij}^n}=\frac{(-1)^n n!t_{ij}}{|\vec{r}_i-\vec{r}_j|^{n+1}}
\f]

