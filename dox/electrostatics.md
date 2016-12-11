Electrostatics                                                 {#electrostatics}
==============

In classical electrostatics, given a field of \f$N\f$ point charges, such that
the \f$j\f$-th one located at \f$\vec{r}_j\f$ and has charge \f$q_j\f$, the
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

