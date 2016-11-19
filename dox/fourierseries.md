Fourier Series                                                        {#fourier}
==============

Torsion and improper torsion angles are typically modled as Fourier series.  In
general a Fourier series has the form:
\f[
E(\theta)=\sum_{n=0}^\infty \frac{1}{2}a_n 
\cos\left(\frac{2\pi n\theta}{P}-\gamma\right),
\f]
(or a form in which cosine is replaced with sine and the values of the
 parameters change to reflect that cosine is sine shifted by 90 degrees).  Here
\f$P\f$ is the periodicity of the function, \f$a_n\f$ is the height of a
particular wave's contribution, and \f$\gamma\f$ is a phase shift.  In modeling
torsions, the angle is periodic every 360 degrees (\f$P=2\pi\f$) and we want the
energy to be 0 at the minimum (i.e. not -1.0), hence we usually add 1 to the
cosine term, the result is:
\f[
E(\theta)=\sum_{n=0}^\infty \frac{1}{2}a_n\left[1+\cos\left(n\theta-\gamma\right)\right],
\f]
Since most atoms involved in force fields make at most four bonds, the
entire potential should be modelable using only the terms up to \f$n=3\f$.  It
is customary to use \f$\gamma\f$ so that a minimum occurs at 180 degrees, which
means for even \f$n\f$, \f$\gamma=2\pi\f$ and for odd \f$n\f$ \f$\gamma=0\f$.

Aside from that the first and second derivatives are easily seen to be:
\f[
\frac{\partial E(\theta)}{\partial \theta}=-\frac{1}{2}\sum_{n=0}^\infty 
n a_n\sin(n\theta-\gamma)
\f]
and:
\f[
\frac{\partial^2 E(\theta)}{\partial^2 \theta}=-\frac{1}{2}\sum_{n=0}^\infty 
n^2 a_n\cos(n\theta-\gamma).
\f]

This pattern continues and we see that the \f$d\f$-th order derivative is:
\f[
\frac{\partial E(\theta)}{\partial \theta}=(-1)^{(d+1)/2}\frac{1}{2}\sum_{n=0}^\infty 
n^d a_n\sin(n\theta-\gamma)
\f]
if \f$d\f$ is odd, and:
\f[
\frac{\partial^2 E(\theta)}{\partial^2 \theta}=(-1)^{d/2}\frac{1}{2}\sum_{n=0}^\infty 
n^d a_n\cos(n\theta-\gamma)
\f]
if \f$d\f$ is even.


## Conventions, i.e. What to Do with the Half?

The Fourier series will ultimately be used for both the torsion and improper
torsion angles.  Unfortunately, the customary half is not present in the form of
the improper torsion.  Furthermore, Tinker style parameter files typically bake
the half into the torsion parameter already anyways.  Unlike the harmonic
oscillator where there is actually an advantage to having the half around (it
simplifies the derivatives), the half of the Fourier series does nothing for us.
For these reasons I have decied to follow Tinker's lead and bake the half into
the torsion amplitude parameter.

## Implementation Details

Keeping with the idea of making the interface as simple as possible, the
FourierSeries class only takes the amplitudes, angles, and periodicities.  This
is because the derivative only depends on these quantities.  The angle is
passed in a somewhat special way, namely as \f$n\theta-\gamma\f$, this is
because this is always the argument to the trig function, regardless of the
derivative order.  Unlike the harmonic oscillator where we can forget about the
details of the shift, for the Fourier series we have to know the value of \f$n\f$
for the derivatives.  Finally, we ask that for torsions with multiple terms, 
i.e. for torsions that have multiple, and \f$n,\gamma,a_n\f$ terms, act as if
each set of parameters is a new torsion (call the function once for each set of
parameters).

