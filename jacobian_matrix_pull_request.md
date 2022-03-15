## Add jacobian_matrix.jl

Hi Carleton and @jacobian contributor. 

Thanks for this incredibly useful and versalite software.

@djturizo and myself are submitting a few proposed fixes and improvements to calc_basic_jacobian_matrix() feature. Most importantly, we have found a bug in the current implementation.

First and foremost, if the admittance matrix is symmetric, the dp/dtheta block of the Power Flow Jacobian should be symmetric and positive definite. This holds for many of the cases but unforunately the current implementation misses this mark a little bit on a few like case14.m.

 also a few new features which I have encapsulated through the file jacobian_matrix.jl. 

Included in this is

- A proposed new struct: JacobianMatrix which is similar to AdmittanceMatrix, which 

In particular, myself and djturizo identified a 