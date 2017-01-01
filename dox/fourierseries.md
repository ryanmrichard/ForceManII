Fourier Series                                                        {#fourier}
==============

Torsion and improper torsion angles are typically modled as Fourier series.  In
general a Fourier series has the form:
\f[
E(\theta)=\sum_{n=0}^\infty \frac{1}{2}a_n 
\cos\left(\frac{2\pi n\theta}{P}-\gamma\right),
\f]
[or a form in which cosine is replaced with sine and the values of the
 parameters change to reflect that \f$\cos(\theta)=\sin(\theta+\pi/2)\f$].  Here
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
Our implementation follows Tinker's convention of absorbing the one-half into
the amplitude.

The first and second derivatives are easily seen to be:
\f[
\frac{\partial E(\theta)}{\partial \theta}=-\sum_{n=0}^\infty
n a_n\sin(n\theta-\gamma)
\f]
and:
\f[
\frac{\partial^2 E(\theta)}{\partial^2 \theta}=-\sum_{n=0}^\infty
n^2 a_n\cos(n\theta-\gamma).
\f]

This pattern continues and we see that the \f$d\f$-th order derivative is:
\f[
\frac{\partial E(\theta)}{\partial \theta}=(-1)^{(d+1)/2}\sum_{n=0}^\infty
n^d a_n\sin(n\theta-\gamma)
\f]
if \f$d\f$ is odd, and:
\f[
\frac{\partial^2 E(\theta)}{\partial^2 \theta}=(-1)^{d/2}\sum_{n=0}^\infty
n^d a_n\cos(n\theta-\gamma)
\f]
if \f$d\f$ is even.

