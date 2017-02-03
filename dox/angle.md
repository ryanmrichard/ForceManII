Angle Class                                                             {#angle}
===========

The Angle class is responsible for normal angle bending terms (i.e. those among
three bonded atoms).  This page describes the mathematical theory behind the
guts of the Angle class.

Mathematical Background
-----------------------

### 0-th Derivative

Assuming the vectors \f$\vec{r_{12}}\f$ and
\f$\vec{r_{32}}\f$ are
vectors parallel to the bond between atoms 1 and 2 and the bond between atoms
2 and 3 respectively the cosine of the angle between these vectors, 
\f$\cos\theta\f$ is given by:

\f[
\cos\theta=\frac{\vec{r}_{12}\cdot\vec{r}_{32}}
                {r_{12} r_{32}},
\f]
where \f$r\f$ denotes the magnitude of \f$\vec{r}\f$.  We can also express this
in terms of the \f$\sin\theta\f$:
\f[
\sin\theta=\frac{||\vec{r}_{12}\times\vec{r}_{32}||}
                {r_{12} r_{32}},
\f]
Consequentially, we can also use the tangent of the angle:
\f[
\tan\theta=\frac{||\vec{r}_{12}\times\vec{r}_{32}||}
                {\vec{r}_{12}\cdot\vec{r}_{32}}.
\f]

The use of the tangent of the angle is preferred numerically as in most coding
languages, there exists a function (usually called `atan2`) that preserves the
sign of the resulting angle and does not suffer from numerical issues at any
angle value (90 degrees, which is not in arctangent's range, is handled
correctly by deffering to cosine at that point).  In particular this means we
avoid the nefarious planar angle problem of most codes!!!!

### Gradient of the Angle

Noting that (as well as defining \f$x=\vec{r}_{12}\times\vec{r}_{32}\f$ and
\f$y=\vec{r}_{12}\cdot\vec{r}_{32}\f$):
\f[
d\theta=d\arctan\left(\frac{x}{y}\right)=\frac{y}{y^2+x^2}dx-\frac{x}{y^2+x^2}dy,
\f]
means we we need two sets of derivatives to obtain
\f$\frac{\partial \theta}{\partial \vec{r}_i}\f$.  The easy one is:
\f{eqnarray*}{
\frac{\partial y}{\partial \vec{r}_1}=&\vec{r}_{32}\\
\frac{\partial y}{\partial \vec{r}_2}=&-\vec{r}_{12}-\vec{r}_{32}\\
\frac{\partial y}{\partial \vec{r}_3}=&\vec{r}_{12}
\f}

Using the chain-rule on the cross product term we get:
\f[
\frac{\partial x}{\partial \vec{r}_i}=
\frac{\partial ||\vec{r}_{12}\times\vec{r}_{32}||}{\partial \vec{r}_i}=
\frac{1}{2||\vec{r}_{12}\times\vec{r}_{32}||}
\frac{\partial \left(\vec{r}_{12}\times\vec{r}_{32}\right)^2}
{\partial \vec{r}_i}=
\frac{\left(\vec{r}_{12}\times\vec{r}_{32}\right)}
     {||\vec{r}_{12}\times\vec{r}_{32}||}\cdot
\frac{\partial\vec{r}_{12}\times\vec{r}_{32}}{\partial \vec{r}_i}
\f]

The cross product can be thought of as a function that takes two vectors and
returns a third vector.  Consequentially, the Jacobian of the cross product is:
\f[
\frac{\partial \vec{r}_{12}\times\vec{r}_{32}}{\vec{r}_{i}}=
\frac{\partial }{\vec{r}_{i}}\begin{bmatrix}
\vec{r}_{12y}\vec{r}_{32z}-\vec{r}_{12z}\vec{r}_{32y}\\
\vec{r}_{12z}\vec{r}_{32x}-\vec{r}_{12x}\vec{r}_{32z}\\
\vec{r}_{12x}\vec{r}_{32y}-\vec{r}_{12y}\vec{r}_{32x}
\end{bmatrix}=
\begin{bmatrix}
\vec{r}_{32z}\frac{\partial \vec{r}_{12y}}{\partial \vec{r}_i}+
\vec{r}_{12y}\frac{\partial \vec{r}_{32z}}{\partial \vec{r}_i}-
\vec{r}_{32y}\frac{\partial \vec{r}_{12z}}{\partial \vec{r}_i}-
\vec{r}_{12z}\frac{\partial \vec{r}_{32y}}{\partial \vec{r}_i}\\
\vec{r}_{32x}\frac{\partial \vec{r}_{12z}}{\partial \vec{r}_i}+
\vec{r}_{12z}\frac{\partial \vec{r}_{32x}}{\partial \vec{r}_i}-
\vec{r}_{32z}\frac{\partial \vec{r}_{12x}}{\partial \vec{r}_i}-
\vec{r}_{12x}\frac{\partial \vec{r}_{32z}}{\partial \vec{r}_i}\\
\vec{r}_{32y}\frac{\partial \vec{r}_{12x}}{\partial \vec{r}_i}+
\vec{r}_{12x}\frac{\partial \vec{r}_{32y}}{\partial \vec{r}_i}-
\vec{r}_{32x}\frac{\partial \vec{r}_{12y}}{\partial \vec{r}_i}-
\vec{r}_{12y}\frac{\partial \vec{r}_{32x}}{\partial \vec{r}_i}&
\end{bmatrix}
\f]

This cleans up nicely for our three vectors of interest:
\f{eqnarray*}{
\frac{\partial \vec{r}_{12}\times\vec{r}_{32}}{\vec{r}_{1}}=&
\begin{bmatrix}
0            &  \vec{r}_{32z} & -\vec{r}_{32y}\\
-\vec{r}_{32z}&  0             &  \vec{r}_{32x}\\
\vec{r}_{32y}&  -\vec{r}_{32x} & 0
\end{bmatrix}\\
\frac{\partial \vec{r}_{12}\times\vec{r}_{32}}{\vec{r}_{2}}=&
\begin{bmatrix}
0                          &\vec{r}_{12z}-\vec{r}_{32z}&\vec{r}_{32y}-\vec{r}_{12y}\\
\vec{r}_{32z}-\vec{r}_{12z}&0                          &\vec{r}_{12x}-\vec{r}_{32x}\\
\vec{r}_{12y}-\vec{r}_{32y}&\vec{r}_{32x}-\vec{r}_{12x}&0
\end{bmatrix}\\
\frac{\partial \vec{r}_{12}\times\vec{r}_{32}}{\vec{r}_{3}}=&
\begin{bmatrix}
0            &  -\vec{r}_{12z} & \vec{r}_{12y}\\
\vec{r}_{12z}&  0             &  -\vec{r}_{12x}\\
-\vec{r}_{12y}&  \vec{r}_{12x} & 0
\end{bmatrix}\\
\f}

Defining \f$\vec{n}=\vec{r}_{12}\times\vec{r}_{32}\f$ this gives:
\f{eqnarray*}{
\frac{\partial x}{\partial \vec{r}_1}=&
\frac{\vec{n}}
     {||\vec{r}_{12}\times\vec{r}_{32}||}\cdot\begin{bmatrix}
     0            &  \vec{r}_{32z} & -\vec{r}_{32y}\\
     -\vec{r}_{32z}&  0             &  \vec{r}_{32x}\\
     \vec{r}_{32y}&  -\vec{r}_{32x} & 0
     \end{bmatrix}=
\frac{1}{||\vec{r}_{12}\times\vec{r}_{32}||}\begin{bmatrix}
\vec{n}_z\vec{r}_{32y}-\vec{n}_y\vec{r}_{32z}\\
\vec{n}_x\vec{r}_{32z}-\vec{n}_z\vec{r}_{32x}\\
\vec{n}_y\vec{r}_{32x}-\vec{n}_x\vec{r}_{32y}\end{bmatrix}\equiv
\frac{\vec{A}}{||\vec{r}_{12}\times\vec{r}_{32}||}\\
\frac{\partial x}{\partial \vec{r}_2}=&
\frac{\vec{n}}
     {||\vec{r}_{12}\times\vec{r}_{32}||}\cdot\begin{bmatrix}
     0            &  -\vec{r}_{31z} & \vec{r}_{31y}\\
     \vec{r}_{31z}&  0             &  -\vec{r}_{31x}\\
     -\vec{r}_{31y}&  \vec{r}_{31x} & 0
     \end{bmatrix}=
\frac{1}{||\vec{r}_{12}\times\vec{r}_{32}||}\begin{bmatrix}
\vec{n}_y\vec{r}_{31z}-\vec{n}_z\vec{r}_{31y}\\
\vec{n}_z\vec{r}_{31x}-\vec{n}_x\vec{r}_{31z}\\
\vec{n}_x\vec{r}_{31y}-\vec{n}_y\vec{r}_{31x}\end{bmatrix}\equiv
\frac{\vec{B}}{||\vec{r}_{12}\times\vec{r}_{32}||}\\
\frac{\partial x}{\partial \vec{r}_3}=&
\frac{\vec{n}}
     {||\vec{r}_{12}\times\vec{r}_{32}||}\cdot\begin{bmatrix}
     0            &  -\vec{r}_{32z} & \vec{r}_{32y}\\
     \vec{r}_{32z}&  0             &  -\vec{r}_{32x}\\
     -\vec{r}_{32y}&  \vec{r}_{32x} & 0
     \end{bmatrix}=
\frac{1}{||\vec{r}_{12}\times\vec{r}_{32}||}\begin{bmatrix}
\vec{n}_y\vec{r}_{32z}-\vec{n}_z\vec{r}_{32y}\\
\vec{n}_z\vec{r}_{32x}-\vec{n}_x\vec{r}_{32z}\\
\vec{n}_x\vec{r}_{32y}-\vec{n}_y\vec{r}_{32x}\end{bmatrix}\equiv
\frac{\vec{C}}{||\vec{r}_{12}\times\vec{r}_{32}||}
\f}

Upon substitution:
\f[
\bigtriangledown\theta=
\frac{1}{\left(||\vec{r}_{12}\times\vec{r}_{32}||\right)^2+
\left(\vec{r}_{12}\cdot\vec{r}_{32}\right)^2}
\begin{bmatrix}
\frac{\vec{A}}{\tan\theta}-\vec{r}_{32}||\vec{r}_{12}\times\vec{r}_{32}||&
\frac{\vec{B}}{\tan\theta}+\left(\vec{r}_{12}+\vec{r}_{32}\right)||\vec{r}_{12}\times\vec{r}_{32}||&
\frac{\vec{C}}{\tan\theta}-\vec{r}_{12}||\vec{r}_{12}\times\vec{r}_{32}||
\end{bmatrix}
\f]

