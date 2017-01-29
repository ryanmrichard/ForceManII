Distance Class                                                       {#distance}
==============

The Distance class is responsible for storing internal coordinates that only
depend on the distance between two points in space.  This includes:
- Bond terms
- Urey-Bradley terms (angle-bond cross terms)
- Charge-charge electrostatic terms
- Dispersion terms
- Exchange terms
 
## Mathematical Description

Given a pair of atoms, the distance \f$r_{12}\f$ between them is:
\f[
   r_{12}=|\mathbf{r_1}-\mathbf{r_2}|= 
           \sqrt{\sum_{i=1}^3\left(q_{1i}-q_{2i}\right)^2},
\f]
where \f$\mathbf{r_i}\f$ is the position of atom \f$i\f$ and \f$q_{ij}\f$ is
the \f$j\f$-th Cartesian coordinate of that atom.
 
The first derivative is:
\f[
    \frac{\partial r_{12}}{\partial q_{ij}}=\frac{1}{2r_{12}}
       \frac{\partial \sum_{k=1}^3\left(q_{1k}-q_{2k}\right)^2}{\partial q_{ij}}
   =\frac{\delta_{i,1}q_{1j}-\delta_{i,2}q_{2j}}{r_{12}}
 \f]

In vector form we have:
\f[
    \frac{\partial r_{12}}{\partial \vec{r}}=
    \left\lbrace\frac{\vec{q}}{r_{12}},\frac{-\vec{q}}{r_{12}}\right\rbrace
\f]

Note that in the case of the angle-bond cross term the internal coordinate that
is usually used is the distance between the ends of the angle.  This does not
depend on the vertex atom, although the parameters do.

The second derivative is:
\f[
    \frac{\partial^2 r_{12}}{\partial q_{ij}\partial q_{kl}}
   =\frac{\partial}{\partial q_{kl}}\frac{\delta_{i,1}q_{1j}-\delta_{i,2}q_{2j}}{r_{12}}
   =\frac{\delta_{l,j}\left(\delta_{i,1}\delta_{k,1}-\delta_{i,2}\delta_{k,2}\right)}{r_{12}}-
    \frac{\delta_{i,1}q_{1j}-\delta_{i,2}q_{2j}}{r_{12}^2}\frac{\partial r_{12}}{\partial q_{kl}}
    =\frac{\delta_{l,j}\left(\delta_{i,1}\delta_{k,1}-\delta_{i,2}\delta_{k,2}\right)}{r_{12}}-
     \frac{\left(\delta_{i,1}q_{1j}-\delta_{i,2}q_{2j}\right)
           \left(\delta_{k,1}q_{1l}-\delta_{k,2}q_{2l}\right)}{r_{12}^3}
 \f]
