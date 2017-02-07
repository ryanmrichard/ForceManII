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
\frac{\partial y}{\partial \vec{r}_1}=&{\vec{r}_{32}}{}^{T}\\
\frac{\partial y}{\partial \vec{r}_2}=&-\left(\vec{r}_{12}+\vec{r}_{32}\right)^T\\
\frac{\partial y}{\partial \vec{r}_3}=&{\vec{r}_{12}}{}^T
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

The derivative distributes over a cross product:
\f[
\frac{\partial\vec{r}_{12}\times\vec{r}_{32}}{\partial \vec{r}_i}=
\frac{\partial\vec{r}_{12}}{\partial \vec{r}_i}\times\vec{r}_{32}+
\vec{r}_{12}\times\frac{\partial\vec{r}_{32}}{\partial \vec{r}_i}
\f]
Noting the identity:
\f[
\vec{a}\cdot(\vec{b}\times\vec{c})=
\vec{b}\cdot(\vec{c}\times\vec{a})=
\vec{c}\cdot(\vec{a}\times\vec{b})
\f]
and defining \f$\vec{n}=\vec{r}_{12}\times\vec{r}_{32}\f$ we now have:
\f[
\vec{n}\cdot\frac{\partial\vec{r}_{12}\times\vec{r}_{32}}{\partial \vec{r}_i}=
\frac{\partial\vec{r}_{12}}{\partial \vec{r}_i}\cdot
\left(\vec{r}_{32}\times\vec{n}\right)+
\frac{\partial\vec{r}_{32}}{\partial \vec{r}_i}\cdot
\left(\vec{n}\times\vec{r}_{12}\right)
\f]
For our three vectors this becomes:
\f{align*}{
\vec{n}\cdot\frac{\partial\vec{r}_{12}\times\vec{r}_{32}}{\partial \vec{r}_1}=&
\vec{r}_{32}\times\vec{n}\\
\vec{n}\cdot\frac{\partial\vec{r}_{12}\times\vec{r}_{32}}{\partial \vec{r}_2}=&
-\vec{r}_{32}\times\vec{n}-\vec{n}\times\vec{r}_{12}=
\vec{n}\times\left(\vec{r}_{32}-\vec{r}_{12}\right)=\vec{n}\times\vec{r}_{31}\\
\vec{n}\cdot\frac{\partial\vec{r}_{12}\times\vec{r}_{32}}{\partial \vec{r}_3}=&
\vec{n}\times\vec{r}_{12}\\
\f}

Upon substitution:
\f[
\bigtriangledown\theta=
\frac{1}{\left(||\vec{r}_{12}\times\vec{r}_{32}||\right)^2+
\left(\vec{r}_{12}\cdot\vec{r}_{32}\right)^2}
\begin{bmatrix}
\frac{\vec{r}_{32}\times\vec{n}}{\tan\theta}-\vec{r}_{32}{}^T||\vec{r}_{12}\times\vec{r}_{32}||&
\frac{\vec{n}\times\vec{r}_{31}}{\tan\theta}+\left(\vec{r}_{12}+\vec{r}_{32}\right)^T||\vec{r}_{12}\times\vec{r}_{32}||&
\frac{\vec{n}\times\vec{r}_{12}}{\tan\theta}-\vec{r}_{12}{}^T||\vec{r}_{12}\times\vec{r}_{32}||
\end{bmatrix}
\f]

