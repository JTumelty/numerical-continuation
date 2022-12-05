# numerical-continuation

MATLAB code to perform numerical continuation on the Swift--Hohenberg equation: 
u_t = ru + (1+d_xx)^2u + nu u^2-u^3, so that solution profiles for steady
state u are found as the parameter r is varied.

To run.
1. create a new directory.
2. initialise a file all_param.mat, that Contains a structure all_param,
    containing a structure ic with fields: r, u, nu, which correspond to the initial
    values of these variables/constants. 
3. Within the current working directory, run contSH_main()

It is possible to specify other parameters (see examples/standard_continuation/all_param.mat and
the file init_param.m for details and information about the default values.)

Notable features include:
- The ability to compute the stability of steady states (set all_param.sp.find_eval = 1)
- The ability to continue bifurcation points in (r,nu) parameter space, given
    initial conditions for r, nu, u and v, where v is the marginal eigenvector
    at the bifurcation.
- The ability to track Maxwell Points, where the free energy is 0, in (r,nu) space. (see examples/mp/)


