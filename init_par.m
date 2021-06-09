function [kp,ic, fp, sd, cp,sp,opt] = init_par()
%INIT_PAR finds the required initial conditions given to the system and
%sets defaults for parameters that are needed but not stated.
%
%INPUTS:Required: 
%           Initial conditions (ic): 
%               Structure ic in file ic.mat in working directory containing
%                   ic.r    Starting value of parameter r
%                   ic.nu   Chosen nu to use in the system
%                   ic.u    Initial profile for a solution to the system
%                           with parameter ic.r
%              
%                   
%
%       Optional:
%           Fixed parameters (fp):
%               Structure fp in file fp.mat in working directory
%                   fp.N    (default 512)   Number of Fourier modes used
%                   fp.Lx   (default 16)    Number of wavelengths in the domain
%                   fp.k                    Wavenumbers in discretised spectral space
%                   fp.x                    Points in discretised real space
%                   fp.w    (default N)     Weight in inner product x =[u,r]
%                                           (x,x) = (u,u) + w(r,r)
%
%           Krylov parameters (kp):
%               Structure kp in file kp.mat in working directory
%                   kp.tolk (default 1e-3)  Relative tolerance in Newton
%                                           method
%                   kp.kmeth (default 0)    Krylov method to solve Jx = F:
%                                                   0 for GMRES
%                                                   1 for BICGStab
%                                                   2 for IDR(s)
%                       (if kp.kmeth = 2
%                       kp.s    (default 20) Dimension of shadow subspace in
%                                            IDR(s).)
%
%                   kp.maxit (default fp.N) Maximum number of iterations in
%                                           Krylov subspace method
%                   kp.precon (default 0)   Preconditioner used: 0 for L^-1 
%                                                                1 for dt(I-dtL)^-1
%                       (if kp.precon = 1
%                       kp.dt   (default 0.1)dt used in preconditioner
%                                            dt(I-dtL)^-1)
%
%           Continuation parameters (cp):
%               Structure cp in file cp.mat in working directory
%                   cp.ds (default 1e-2)    Initial step size in
%                                           continuation process
%                   cp.maxitn (default 30)  Maximum number of Newton
%                                           iterations
%                   cp.toln (default 1e-9)  Tolerance for convergence of
%                                           Newton iteration
%                   cp.ts   (default 10)    Store results every cp.ts steps
%                                           in continuation process
%                   cp.Ns   (default 1000)  Number of steps taken in
%                                           continuation process
%                   cp.vds  (default 0)    Vary step size: 0 for no
%                                                           1 for yes
%                       (if cp.vds = 1
%                       cp.minds (default = 1e-6) Minimum allowed step size
%                       cp.maxds (default = 0.01) Maximum allowed step size
%                       cp.mincost (optional)     Minimum cost of iterations
%                                                 leading to increase the stepsize
%                       cp.maxcost (optional)     Maximum cost of iterations
%                                                 leading to decrease the
%                    	                          stepsize)
%		    cp.mp (defualt 0)	    Follow the Maxwell point: 0 for no
%								      1 for yes
%	            cp.bt (default 0)       Follow a bifurcation point: 0 for no
%									1 for yes
%
%           Stability parameters (fp):
%               Structure sp in file sp.mat in working directory
%               sp.find_eval (default 0) Find eigenvalues?   0 for no
%                                                               1 for yes
%                   (if sp.find_eval = 1
%                   sp.Smeth (default 0)    Time-step scheme to find
%                                           eigenvalues:    0 ETDRK4
%                                                           1 ETDRK2
%                                                           2 Euler
%                   sp.bps (default 10)   Compute eigenvalues every ...
%                                         steps
%                   sp.Sdt (default 0.01 (ETDRK4),0.001 (ETDRK2), 0.0001 (Euler))
%                                         Time-step used to find
%                                         eigenvalues
%                   sp.K (default 20 (ETDRK4), 30 (ETDRK2), 30 (Euler))
%                                          Krylov subspace dimension used
%                   sp.eigstol (default 1e-8) Tolerance used when computing
%                                             eigenvalues
%                   sp.maxit (default 800)  Maximum number of iterations
%                   sp.evalno (default 8)  Number of leading eigenvalues to
%                                          compute)
%
%        Optional parameters (opt):
%               Structure opt in file opt.mat in working directory
%                   opt.armijo (default 0)  Use Armijo rule? 0 for no
%                                                            1 for yes
%                       (if opt.armijo = 1)
%                       opt.aalpha (default 1e-4) (|Jx -F|< (1 - aalpha2^-l)|F|)
%                       opt.maxsteps (default 100) Maximum number of decreases
%                                           in correction size)
%                   opt.combtol (default 0) Combine absolute and relative
%                                           tolerance in Krylov iteration? 
%                                                            0 for no, 
%                                                            1 for yes
%                       (if opt.combtol = 1)
%                       opt.ta (default 1e-3)   Absolutre tolerance
%                       opt.tr (default 1e-6)    Relative tolerance)
%                   opt.plot_profile        A list of iterations to plot
%                                           the profile u
%                                          
%OUTPUTS:   kp
%           ic
%           fp
%           sd
%           cp
%           sp
%           opt

% Find working directory and make log and error files. 
if isfile(strcat(pwd,'/sd.mat'))
    load(strcat(pwd,'/sd.mat'),'sd');
end

sd.dir = pwd;
if ~isfield(sd,'log')
    sd.log = strcat(pwd,'/o.log');
end
if ~isfield(sd,'dataname')
    sd.dataname     = 'data';
end

% Load the mat files containing the initial conditions and parameters.
try
    load(strcat(sd.dir,'/ic.mat'),'ic');
catch
    error('Error. ic not found. Check that ic.mat containing ic exists.')
end
for par = ["fp","cp","sp","kp","opt"]
    try
        load(strcat(sd.dir,'/',par,'.mat'),par);
    catch
        continue
    end
end

% Checks that the required initial conditions (u,r,nu) are given
    if ~isfield(ic,'u')
        error('Error. u not found in initial conditions')
    elseif ~isfield(ic,'r')
        error('Error.  r not found in initial conditions')
    elseif ~isfield(ic,'nu')
        error('Error. nu not found in initial conditions')
    end
    
% Set defaults for fixed parameters (fp)
    if ~exist('fp','var')
        fp = struct;
    end
    
    if ~isfield(fp,'N')
        fp.N            = 512;
    end
    if ~isfield(fp,'Lx')
        fp.Lx           = 16;
    end
    fp.k            = [0:fp.N/2-1 -fp.N/2:-1]'/fp.Lx;    
    fp.x            = -fp.Lx*pi + (2*pi*fp.Lx)/fp.N*(0:(fp.N-1));
    if ~isfield(fp,'w')
        fp.w            = fp.N;
    end

% Set defaults for Krylov iteration process (kp)
    if ~exist('kp','var')
        kp = struct;
    end    
        if ~isfield(kp,'tolk')
            kp.tolk         = 1e-3;
        end
        if ~isfield(kp,'kmeth') 
            % Choose Krylov method: 0 GMRES, 1 BICGStab, 2 IDR(s)
            kp.kmeth        = 0;
        end
        if kp.kmeth == 2 && ~isfield(kp,'s')
            kp.s           = 20;
        end
        if ~isfield(kp,'maxit')
            kp.maxit         = fp.N;
        end
        if ~isfield(kp,'precon')
            kp.precon        = 0; % Choose precondtioner to use: 0 L^-1, Pdt where Pdt(I-PdtL)^-1 
        end

        if kp.precon == 1 && ~isfield(kp,'dt')
            kp.dt           = 0.1;
        end
        

% Set parameters for the continuation process (cp)
    if ~exist('cp','var')
        cp = struct;
    end	
        if ~isfield(cp,'ds')
            cp.ds           = 0.01;
        end
        if ~isfield(cp,'maxitn')
            cp.maxitn       = 30;
        end
        if ~isfield(cp,'toln')
            cp.toln         = 1e-9;
        elseif cp.toln <= 1e-10
            warning('Tolerance for Newton iteration may be too small for some solutions.')
        end
        if ~isfield(cp,'ts')
            cp.ts           = 10;
        end
        if ~isfield(cp,'Ns')
            cp.Ns           = 1000;
        end
        if ~isfield(cp,'vds')
            cp.vds          = 0;
        end
        if cp.vds == 1
            if ~isfield(cp,'minds')
                cp.minds        = 1e-6;
            end
            if ~isfield(cp,'maxds')
                cp.maxds        = 0.01;
            end
        end
        if ~isfield(cp,'bt')
            cp.bt          = 0;
        end
        if ~isfield(cp,'mp')
            cp.mp          = 0;
        end
% Set parameters for computing stability (sp)
    if ~exist('sp','var')
        sp = struct;
    end
        if ~isfield(sp,'find_eval')
            sp.find_eval = 0;
        end

        if sp.find_eval == 1
            if ~isfield(sp,'bps')
                sp.bps          = 10;
            end
            if ~isfield(sp,'Smeth')
                sp.Smeth        = 0;
            end
            if ~isfield(sp,'Sdt')
                if sp.Smeth == 0
                    sp.Sdt = 1e-2;
                elseif sp.Smeth == 1
                    sp.Sdt = 1e-3;
                elseif sp.Smeth == 2
                    sp.Sdt = 1e-4;
                end
            end

            if ~isfield(sp,'K')
                if sp.Smeth == 0
                    sp.K  = 20;
                elseif sp.Smeth == 1
                    sp.K = 30;
                elseif sp.Smeth == 2
                    sp.K = 30;
                end
            end
            if ~isfield(sp,'eigstol')
                sp.eigstol      = 1e-8;
            end
            if ~isfield(sp,'eigsmaxit')   
                sp.eigsmaxit    = 800;
            end
            if ~isfield(sp,'evalno')
                sp.evalno       = 8;
            end
        end

% Set parameters for further optional processes (opt)
    if ~exist('opt','var')
        opt = struct;
    end
    
        if ~isfield(opt,'armijo')
            opt.armijo = 0;
        end
        if opt.armijo == 1 && ~isfield(opt,'aalpha')
            opt.aalpha = 1e-4;
        end
        if opt.armijo == 1 && ~isfield(opt,'maxsteps')
            opt.maxsteps = -100;
        end

        if ~isfield(opt, 'combtol')
            opt.combtol = 0;
        end
        if opt.combtol == 1 && ~isfield(opt,'tr')
            opt.tr = 1e-3;
        end
        if opt.combtol == 1 && ~isfield(opt,'ta')
            opt.ta = 1e-6;
        end
        
        if ~isfield(opt,'plot_profile')
            opt.plot_profile = 0;
        end
        
end

