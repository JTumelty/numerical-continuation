# numerical-continuation

MATLAB code to perform numerical continuation on the Swift--Hohenberg equation: u_t = ru + (1+d_xx)^2u + nu u^2-u^3. 
To run, create a new directory and initialise with files fp.mat, kp.mat, cp.mat and ic.mat that contain a structure of the same name (see test case in the folder gmres_pat_Lx32). The first three structures may be empty, but ic in ic.mat must contain fields: r (the initial value of the bifurcation parameter), u (an initial guess of the steady solution at r) and nu (a fixed value of the parameter in the Swift--Hohenberg equation). Then run contSH_main.m.
