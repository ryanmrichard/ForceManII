Angle Class                                                             {#angle}
===========

The angle class is responsible for normal angle bending terms (i.e. those among
three bonded atoms).  Assuming the vectors \f$\vec{r_{12}}\f$ and 
\f$\vec{r_{32}}\f$ are
vectors parallel to the bond between atoms 1 and 2 and the bond between atoms
2 and 3 respectively the cosine of the angle between these vectors, 
\f$\cos\theta\f$ is given by:

\f[
\cos\theta=\frac{\vec{r_{12}}\cdot\vec{r_{32}}}
                {r_{12} r_{32}},
\f]
where \f$r\f$ denotes the magnitude of \f$\vec{r}\f$.  We can also express this
in terms of the \f$\sin\theta\f$:
\f[
\sin\theta=\frac{\left(\vec{r_{12}}\times\vec{r_{32}}\right)\cdot\vec{n}}
                {r_{12} r_{32} n},
\f]
where \f$\frac{\vec{n}}{n}\f$ is a unit vector perpendicular to both 
\f$\vec{r_{12}}\f$ and \f$\vec{r_{32}}\f$.  Consequentially, we can also use
the tangent of the angle:
\f[
\tan\theta=\frac{\left(\vec{r_{12}}\times\vec{r_{32}}\right)\cdot\vec{n}}
                {\vec{r_{12}}\cdot\vec{r_{32}}n}.
\f]

The use of the tangent of the angle is preferred numerically as in most coding
languages, there exists a function (usually called `atan2`) that preserves the
sign of the resulting angle and does not suffer from numerical issues at any
angle value (90 degrees, which is not in arctangent's range is handled
correctly, as is the nefarious angle of 180 degrees which is outside of
arccosine's range).
