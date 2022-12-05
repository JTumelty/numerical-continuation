function all_res = init_all_res(all_param,all_res)
%INIT_RES initialises the structures all_res containing res, nres, stab_res,
% which hold all results, results for a single continuation step and 
% stability results respectively.
%
%INPUTS:    all_param
%           (optional)
%               all_res: Presence of a second input argument indicates that only
%                        res and stab_res should be initialised and
%                        res.res_no should be incremented
%
%OUTPUTS:   all_res

% Pre-generate the dimension of the result arrays to be the max of
% cp.n_it_to_save, or cp.Ns
size = min(all_param.cp.n_it_to_save,all_param.cp.Ns);

if nargin == 2
    % Store the current res_no
    curr_res_no = all_res.res.res_no;
else
    % Set default number of results stored to 0
    curr_res_no = 0;
    
    % Only initialise nres at the start of the program.
    all_res.nres = init_nres(   all_param.ic,...
                                all_param.cp.ds, ...
                                all_param.cp.maxitn, ...
                                all_param.fp.N,...
                                all_param.cp.mp,...
                                all_param.cp.bt);
end
all_res.res = init_res(     all_param.fp.N, ...
                            all_param.cp, ...
                            size,...
                            curr_res_no);
                        
% Only creates stab_res when we compute stability
if all_param.sp.find_eval == 1
    all_res.stab_res = init_stab_res(   all_param.fp.N,...
                                        size);
end
end
%-------------------------------------------------------------------------%
%% Initialise results
function res = init_res(N,cp,size,curr_res_no)
    %INIT_RES initialises the structure res, which contains details about
    %the current continuation step.
    %INPUTS:    N       dimension of the system
    %           maxitn  Maximum number of Newton iterations
    %           size    Maximum number of results to store before saving
    %          
    %OUTPUTS:   res
    
    % Solution profile and bifurcation parameter
    res.u           = zeros(N,size);
    res.r           = zeros(1,size);
    
    % Save r and E at every continuation step
    res.rE          = zeros(2,size*cp.ts);
    
    % Extra variables if performing non-standard continuation
    if cp.mp == 1
        res.nu      = zeros(1,size);
    elseif cp.bt == 1
        res.v       = zeros(N,size);
        res.nu      = zeros(1,size);
    end
    
    % Global variables to characterise the result
    res.E           = zeros(1,size);
    res.F           = zeros(1,size);
    res.FL          = zeros(1,size);

    % Details about the arc-length along the branch
    res.s           = zeros(1,size);
    res.ds_all      = zeros(1,size);
    
    % Properties indicating the performance of each continuation step
    res.f           = zeros(1,size);
    res.fall        = zeros(cp.maxitn+1, size);
    res.cost        = zeros(1,size);
    
    % Counter for index to update in res
    res.sc          = 1;
    res.rE_sc       = 1;
    
    % Maximum number of results to store.
    res.maxsc       = size;
    % Increment number of results previously stored
    res.res_no      = curr_res_no + 1;
end

%% Initialise current iterate results
function nres = init_nres(ic,ds,maxitn,N,mp,bt)
    %INIT_NRES initialises the structure nres, which contains details about
    %the current continuation step.
    %INPUTS:    ic    initial conditions
    %           ds    initial step size
    %           N     dimension of the system
    %
    %OUTPUTS:   nres
    
    % Initialise the continuation step with the initial conditions
    nres.un         = ic.u;
    nres.rn         = ic.r;
    nres.nun        = ic.nu;
    
    % Initialise the branch length and starting step size
    nres.sn         = 0;
    nres.ds         = ds;
    
    % Initial global variables
    nres.E           = 0;
    nres.F           = 0;
    nres.FL          = 0;
    
    % Store the previous two iterates for prediction steps
    nres.uprev      = [zeros(N,1) ic.u];
    nres.rprev      = [0 ic.r];
    nres.sprev      = [0 0];

    % Iteration number
    nres.its        = 1;
    
    % Initialise variables about Newton iterations
    nres.itn        = 0;
    nres.restart    = 0;
    nres.f          = 0;
    nres.cost       = 0;
    nres.fall       = zeros(maxitn+1,1);
    nres.converged  = 0;
    nres.flag       = 0;
    
    % Add extra initial conditions if performing non-standard continuation
    if mp == 1
        nres.nuprev = [0 ic.nu];
    elseif bt == 1
        nres.vn     = ic.v;
        nres.vprev  = [zeros(N,1) ic.v];
        nres.nuprev = [0 ic.nu];
    end
    
end

%% Initialise stability results
function stab_res = init_stab_res(N,size)
    
    % Initialise array to store (r,E) pairs for points when the eigenvalues
    % and eigenvectors are stored
    stab_res.evalrE         = zeros(2,size);
    
    % Arrays for eigenvectors and eigenvalues
    stab_res.evec           = zeros(sp.evalno*(N),size);
    stab_res.eval           = zeros(sp.evalno,size);

    % Counter for the index to store the eigenvector/eigenvalue results
    stab_res.bp_sc          = 1;
end