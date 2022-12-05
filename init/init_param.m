function all_param = init_param()
%INIT_PAR finds the required initial conditions given to the system and
%sets defaults for parameters that are needed but not stated. Throws an
%error if unexpected format.
%
%EXPECTATION:
%       There should be the file all_param.mat in the current working
%       directory, which containes the structure:
%           all_param                   All parameters
%               with substructures:
%                   Required:   ic      Initial conditions
%                   Optional:   fp      Fixed parameters
%                               cp      Continuation parameters
%                               kp      Krylov parameters
%                               sp      Stability parameters
%                               op      Optional parameters                              
%                               sd      Save directory
%OUTPUTS:
%       all_param: 

    % Check that the correct form exists
    if check_all_param_file_exists()
        all_param = check_all_param_format();
    end

    % Validates that initial conditions were specified and they are valid.
    check_ic_exists(all_param);
    check_ic_valid(all_param.ic)

    % Sets default values for fp, cp, kp,sp, op and sd
    fp = get_param(all_param, 'fp');
    all_param.fp = set_fp(fp, size(all_param.ic.u));
    
    cp = get_param(all_param, 'cp');
    all_param.cp = set_cp(cp);
    
    kp = get_param(all_param, 'kp');
    all_param.kp = set_kp(kp, all_param.fp.N);
    
    sp = get_param(all_param, 'sp');
    all_param.sp = set_sp(sp);
    
    op = get_param(all_param, 'op');
    all_param.op = set_op(op);
    
    sd = get_param(all_param, 'sd');
    all_param.sd = set_sd(sd);
    
    if all_param.cp.bt == 1
        check_ic_valid_bt(all_param.ic);
    end
end

%% Validation of all_param.mat
function successful = check_all_param_file_exists()
    % Confirms that there is a file all_param.mat in the current working
    % directory.
    if ~isfile(strcat(pwd,'/all_param.mat'))
        error('Error. Expected file all_param.mat in the current working directory.')
    end
    successful = true;    
end

function all_param = check_all_param_format()
    %Checks that the mat file all_param.mat is valid and contains a
    %structure all_param.
    try
        all_param = whos('-file','all_param.mat','all_param');
    catch
        error('Error. all_param.mat is not a valid .mat file.')
    end

    if size(all_param,1) == 0
       error('Error. all_param.mat does not contain the structure all_param');
    else
        load('all_param.mat','all_param');
    end
end

%% Validation and setup of ic (Initial conditions) in all_param
% Initial conditions (ic): 
%   Required structure in all_param with:
%       ic.r    Starting value of bifurcation parameter r
%       ic.nu   Initial values of parameters in the system 
%       ic.u    Initial profile for a solution to the system
%               with parameter ic.r
%       
%   If performing bifurcation tracking
%       ic.v    Initial profile for the eigenvecotr of the bifurcation
%               point

function check_ic_exists(all_param)
    % Confirms that all_param contains the required structure ic.
    if ~isfield(all_param,'ic')
        error('Error. all_param does not contain a structure ic')
    end
end

function check_ic_valid(ic)
% Checks that the structure ic contains valid initial conditions (u,r,nu)

    if ~isfield(ic,'u')
        error('Error. u not found in initial conditions')
    elseif ~isa(ic.u,'numeric')
        error('Error. u does not have the correct type. Should be numeric.')
    end

    if ~isfield(ic,'r')
        error('Error.  r not found in initial conditions')
    elseif ~isa(ic.r,'numeric')
        error('Error. r does not have the correct type. Should be numeric.')
    end
    
    if ~isfield(ic,'nu')
        error('Error. nu not found in initial conditions')
    elseif ~isa(ic.nu,'numeric')
        error('Error. nu does not have the correct type. Should be numeric.')
    end
end

function check_ic_valid_bt(ic)
% Checks that the structure ic contains an initial condition for an
% eigenvector v needed for bifurcation tracking

    if ~isfield(ic,'v')
        error('Error. For bifurcation tracking, also need an initial condition v for the eigenvector')
    elseif ~isa(ic.v,'numeric')
        error('Error. v does not have the correct type. Should be numeric.')
    end
end

%% Utility functions
function param = set_defaults(param, default_param)
    %SET_DEFAULTS sets the values of fields in param to the corresponding
    %field values in default_param.
    fields = fieldnames(default_param);
    for i = 1:length(fields)
        % Only sets the field if not given
        if ~isfield(param,fields{i})
            param = setfield(param,fields{i},getfield(default_param,fields{i}));
        end
    end
end

function param = get_param(all_param, param_name)
    %GET_SUB_PARAM gets the field with name param_name from all_param if it
    %exists, otherwise returns empty struct.
    if isfield(all_param,param_name)
        param = getfield(all_param, param_name);
    else
        param = struct;
    end  
end

%% Validation and setup of optional fp (fixed parameters) in all_param
% Fixed parameters (fp):
%   Optional structure fp in all_param with:
%       fp.N    (default 512)   Number of Fourier modes used
%       fp.Lx   (default 16)    Number of wavelengths in the domain
%       fp.k                    Wavenumbers in discretised spectral space
%       fp.x                    Points in discretised real space
%       fp.w    (default fp.N)  Weight in inner product x =[u,r]
%                               (x,x) = (u,u) + w(r,r)

function default_fp = get_default_fp()
%Create the structure default_fp that specifies the default values for
%fields in the structure fp.
    default_fp = struct('N',512,...
                        'Lx',16);
end

function fp = set_fp(fp,N)
    % Set defaults for fixed parameters (fp). N is the dimension of ic.u    
    
    % Set defaults
    fp = set_defaults(fp,get_default_fp());
    
    % Raise an error if the dimension of ic.u does not equal fp.N
    if  fp.N ~= N
        error('Error. fp.N should have the same dimension as ic.u')
    end
    
    % Set the default for fp.w to equal fp.N
    if ~isfield(fp,'w')
        fp.w  = fp.N;
    end
    
    % Construct vectors for k and x directly from N and Lx
    fp.k = [0:fp.N/2-1 -fp.N/2:-1]'/fp.Lx;    
    fp.x  = -fp.Lx*pi + (2*pi*fp.Lx)/fp.N*(0:(fp.N-1));

end

%% Validation and setup of optional kp (Krylov parameters) in all_param
% Krylov parameters (kp):
%   Optional structure kp in all_param: with
%       kp.tolk (default 1e-3)  Relative tolerance in Newton
%                               method
%       kp.kmeth (default 0)    Krylov method to solve Jx = F:
%                                       0 for GMRES
%                                       1 for BICGStab
% 
%       kp.maxit (default fp.N) Maximum number of iterations in
%                               Krylov subspace method
%       kp.precon (default 0)   Preconditioner used: 0 for L^-1 
%                                                    1 for dt(I-dtL)^-1
%        if kp.precon = 1
%           kp.dt   (default 0.1)dt used in preconditioner
%                                dt(I-dtL)^-1
%

function default_kp = get_default_kp()
%Create the structure default_kp that specifies the default values for
%fields in the structure kp.
    default_kp = struct('tolk',1e-3,...
                        'kmeth',0,...
                        'precon',0,...
                        'dt',0.1);
end

function kp = set_kp(kp,N)
    % Set defaults for Krylov iteration process (kp)

    % Set defaults
    kp = set_defaults(kp,get_default_kp());

    % Set dependent parameters
    if ~isfield(kp,'maxit')
        kp.maxit = N;
    end
end

%%  Validation and setup of optional cp (Continuation parameters) in all_param
% Continuation parameters (cp):
%   Optional structure cp in all_param with:
%       cp.ds (default 1e-2)    Initial step size in
%                               continuation process
%       cp.maxitn (default 30)  Maximum number of Newton
%                               iterations
%       cp.toln (default 1e-9)  Tolerance for convergence of
%                               Newton iteration
%       cp.ts   (default 10)    Store results every cp.ts steps
%                               in continuation process
%       cp.Ns   (default 1000)  Number of steps taken in
%                               continuation process
%       cp.vds  (default 0)    Vary step size: 0 for no
%                                               1 for yes
%       if cp.vds = 1:
%           cp.minds (default = 1e-6) Minimum allowed step size
%           cp.maxds (default = 0.01) Maximum allowed step size
%           cp.mincost (optional)     Minimum cost of iterations
%                                     leading to increase the stepsize
%           cp.maxcost (optional)     Maximum cost of iterations
%                                     leading to decrease the
%        	                          stepsize)
%           cp.mp (defualt 0)	      Follow the Maxwell point: 0 for no
% 					                  1 for yes
%           cp.bt (default 0)         Follow a bifurcation point: 0 for no
% 						              1 for yes
%

function default_cp = get_default_cp()
%Create the structure default_cp that specifies the default values for
%fields in the structure cp.
    default_cp  = struct(   'ds',0.01,...
                            'maxitn',30,...
                            'ts',10,...
                            'toln',1e-9,...
                            'Ns',1000,...
                            'vds',0,...
                            'minds',1e-6,...
                            'maxds',0.01,...
                            'bt',0,...
                            'mp',0,...
                            'n_it_to_save',100);
end

function cp = set_cp(cp)
    % Set parameters for the continuation process (cp)
    
    % Set defaults
    cp = set_defaults(cp,get_default_cp());
    
    % Provide warning if the tolerance is too small
    if cp.toln <= 1e-10
        warning('Tolerance for Newton iteration may be too small for some solutions.')
    end
    
    % Check continuation options are valid
    if cp.bt == 1 && cp.mp == 1
        error('Error. Cannot perform both bifurcation tracing and continuation of the Maxwell point.') 
    end
end

%% Validation and setup of optional sp (Stability parameters) in all_param
% Stability parameters (sp):
%    Optional structure sp in fall_param with:
%       sp.find_eval (default 0) Find eigenvalues?   0 for no
%                                                    1 for yes
%       
%           if sp.find_eval = 1:
%               sp.Smeth (default 0)    Time-step scheme to find
%                                       eigenvalues:    0 ETDRK4
%                                                       1 ETDRK2
%                                                       2 Euler
%               sp.bps (default 10)   Compute eigenvalues every ...
%                                     steps
%               sp.Sdt (default 0.01 (ETDRK4),0.001 (ETDRK2), 0.0001 (Euler))
%                                     Time-step used to find
%                                     eigenvalues
%               sp.K (default 20 (ETDRK4), 30 (ETDRK2), 30 (Euler))
%                                      Krylov subspace dimension used
%               sp.eigstol (default 1e-8) Tolerance used when computing
%                                         eigenvalues
%               sp.maxit (default 800)  Maximum number of iterations
%               sp.evalno (default 8)  Number of leading eigenvalues to
%                                      compute

function default_sp = get_default_sp()
%Create the structure default_sp that specifies the default values for
%fields in the structure sp.
    default_sp = struct('find_eval',0,...
                        'bps',10,...
                        'Smeth',0,...
                        'eigstol',1e-8,...
                        'eigsmaxit',800,...
                        'evalno',8);
end

function sp = set_sp(sp)
    % Set parameters for computing stability (sp) 
    
   sp = set_defaults(sp,get_default_sp());
   
   %Default values for Sdt and K depend on Smeth, so should be set
   %separately
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
end

%% Validation and setup of optional op (optional parameters) in all_param
%        Optional parameters (op):
%               Structure op in all_param
%
%                   op.combtol (default 0) Combine absolute and relative
%                                           tolerance in Krylov iteration? 
%                                                            0 for no, 
%                                                            1 for yes
%                            (if op.combtol = 1)
%                                 op.ta (default 1e-3)   Absolutre tolerance
%                                 op.tr (default 1e-6)    Relative tolerance)
%
%                   op.plot_profile        A list of iterations to plot
%                                           the profile u
%                   op.n_it_to_save    (default 10000) Save results at max
%                                         every n_it_to_save steps
%                                          

function default_op = get_default_op()
%Create the structure default_op that specifies the default values for
%fields in the structure op.
    default_op = struct('combtol',0,...
                        'tr',1e-3,...
                        'ta',1e-6,...
                        'plot_profile',0);

end

function op = set_op(op)
   % Set parameters for further optional processes (op)
    
    % Set defaults
    op = set_defaults(op,get_default_op());

end

%% Validation and setup of optional sd (save directory) in all_param
%        Save directory (sd):
%               Structure sd in all_param
%                   sd.dir:         The directory to save the results. Default:
%                                   current working directory
%                   sd.dataname:    The filename for the results to be saved in     
%                   sd.log:         Boolean for whether to log the results
%                   sd.log_name:    File name to log the results under

function default_sd = get_default_sd()
%Create the structure default_op that specifies the default values for
%fields in the structure op.
    default_sd = struct('dir',pwd,...
                        'dataname','/data',...
                        'log',1);

end

function sd = set_sd(sd)
   % Set parameters for further save directory (sd)
   
    % Set defaults
    sd = set_defaults(sd,get_default_sd());
    
    % Set default path for log file
    if ~isfield(sd,'log_name')
        sd.log_name = strcat(sd.dir,'/o.log');
    end
    
    
end