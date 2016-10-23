Distance Class                                                       {#distance}
==============

The Distance class is responsible for storing internal coordinates that only
depend on the distance between two points in space.  This includes:
- Bond terms
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
   =\frac{\pm \left(q_{1j}-q_{2j}\right)}{r_{12}}
 \f]
where the plus is for \f$i=1\f$ and the minus for \f$i=2\f$.
