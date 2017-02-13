Harmonic Oscillator                                                        {#HO}
===================
![](ho.svg)


One of the cornerstones of physics' approximations is the harmonic oscillator.
The energy of an \f$N\f$-dimensional harmonic oscillator is given by:
\f[
    E(\vec{q};\vec{k},\vec{q}_0)=\frac{1}{2}\vec{k}\cdot
    \left[\left(\vec{q}-\vec{q}_0\right)\circ\left(\vec{q}-\vec{q}_0\right)\right],
\f]
where \f$\vec{q}\f$ and \f$\vec{q}_0\f$ are vectors of \f$N\f$ coordinates such
that \f$q_i\f$ and \f$[q_0]_i\f$ respectively are the current and equilibrium
values of the \f$i\f$-th coordinate.  \f$\vec{k}\f$ is a
vector of \f$N\f$ force constants such that the \f$q_i\f$ has force constant 
\f$k_i\f$.  \f$\circ\f$ denotes the Hadamard product (element-wise
multiplication), i.e. \f$\left[\vec{q}\circ\vec{q}\right]_{i}=q_i^2\f$.  In
element-wise notation this reduces to the more familiar form:
\f[
    E(\vec{q};\vec{k},\vec{q}_0)=\sum_{i=1}^N\frac{1}{2}k_i
    \left(q_i-\left[q_0\right]_i\right)^2,
\f]

The derivative w.r.t. to \f$q_i\f$ is:
\f[
   \frac{\partial E(\vec{q};\vec{k},\vec{q}_0)}{\partial q_i}=
     k_i\left(q_i-\left[q_0\right]_i\right),
\f]
the Hessian is just \f$\vec{k}\f$ down the diagonal.  All higher derivatives 
are 0.  Doesn't get much easier than that.

## Notes on H.O. as Implemented in ForceManII

- As with all force field quantities we assume atomic units are used throughout
- We explicitly include the 1/2 in our energy formula as it simplifies the
  derivatives
