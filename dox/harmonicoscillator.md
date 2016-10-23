Harmonic Oscillator                                                        {#HO}
===================

One of the cornerstones of physics' approximations is the harmonic oscillator.
The energy of an \f$N\f$-dimensional harmonic oscillator is given by:
\f[
    E(\vec{q};\vec{k})=\frac{1}{2}\vec{k}\cdot(\vec{q}\circ\vec{q}),
\f]
where \f$\vec{q}\f$ is a vector of \f$N\f$ coordinates and \f$\vec{k}\f$ is a 
vector of \f$N\f$ force constants such that the \f$q_i\f$ has force constant 
\f$k_i\f$.  \f$\circ\f$ denotes the Hadamard product (element-wise
multiplication), i.e. \f$\left[\vec{q}\circ\vec{q}\right]_{i}=q_i^2\f$.

The derivative w.r.t. to \f$\vec{q}\f$ is:
\f[
   \left[\frac{\partial E(\vec{q};\vec{k})}{\partial\vec{q}}\right]_i=k_i q_i,
\f]
the Hessian is just \f$\vec{k}\f$ down the diagonal.  All higher derivatives 
are 0.  Doesn't get much easier than that.

## Notes on H.O. as Implemented in ForceManII

- As with all force field quantities we assume atomic units are used throughout
- We explicitly include the 1/2 in our energy formula
- Usually a H.O. measures the energy penelty caused by displacemet from some
  value \f$\vec{q_0}\f$.  To simplify the underlying code we assume that input
  coordinates \f$\vec{q}=\vec{q}^\prime-\vec{q_0}\f$ where \f$\vec{q}^\prime\f$
  is the actual value of the bond, angle, etc.
