Polynomial Potential Terms                                         {#polynomial}
==========================

\note The polynomial term is not yet implemented in ForceManII, but will be at
some point.

Many of the forms of force field potentials are that of polynomials.  That is
they are a linear combination of \f$N\f$ terms of the form:
\f[
E=\sum_{i=1}^Nc_ix^{p_i}
\f]
where \f$c_i\f$ is the weight of the \f$i\f$-th term which has the exponet 
\f$p_i\f$.  In general \f$p_i\f$ can be positive (say \f$p_i=2\f$ for the
harmonic oscillator) or negative (say \f$p_i=-1\f$ for point-charge,
point-charge interactions).  The gradient is:
\f[
\frac{\partial E}{\partial x}=\sum_{i=1}^N c_ip_i x^{p_i-1}
\f]
and the Hessian is:
\f[
\frac{\partial^2 E}{\partial x^2}=\sum_{i=1}^N c_ip_i(p_i-1)x^{p_i-2}
\f]
At the moment the Harmonic Oscillator term uses its own class.
