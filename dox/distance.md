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
the \f$j\f$-th Cartesian coordinate of that atom.  In vector format this is:
\f[
r_{12}=\sqrt{\vec{r}_{1}^2-2\vec{r}_1\cdot\vec{r}_2+\vec{r}_2^2}
\f]
 
The derivative with respect to \f$r_1\f$ is:
\f[
\frac{\partial r_{12}}{\partial \vec{r}_1}=
\frac{1}{2r_{12}}2\vec{r}_1-2\vec{r}_2=\frac{\vec{r}_{12}}{r_{12}}
\f]
and the derivative with respect to \f$r_2\f$ is the negative of this.  Therefore
the gradient is:
\f[
\bigtriangledown r_{12}=\frac{1}{r_{12}}
   \begin{bmatrix}\vec{r}_{12} & -\vec{r}_{12}\end{bmatrix}
\f]

The gradient of the gradient is:
\f[
\bigtriangledown^2 r_{12}=\bigtriangledown (\bigtriangledown r_{12})^T=
\bigtriangledown\frac{1}{r_{12}}
   \begin{bmatrix}\vec{r}_{12}\\-\vec{r}_{12}\end{bmatrix}=
   \frac{-\left(\bigtriangledown r_{12}\right)^T}
        {r_{12}}\cdot\left(\bigtriangledown r_{12}\right)+
   \frac{1}{r_{12}}
        \bigtriangledown\begin{bmatrix}\vec{r}_{12}\\-\vec{r}_{12}\end{bmatrix}=
   \frac{1}{r_{12}}\left(\begin{bmatrix}\mathbf{1}&-\mathbf{1}\\
   -\mathbf{1}&\mathbf{1}\end{bmatrix}-
   \left(\bigtriangledown r_{12}\right)^T\cdot\bigtriangledown r_{12}\right)
\f]

where a superscript \f$T\f$ denotes a transpose and \f$\mathbf{1}\f$ is the 3 by
3 identity matrix.

\note In the case of the angle-bond cross term the internal coordinate that
is usually used is the distance between the ends of the angle.  This does not
depend on the vertex atom, although the parameters do.
