function  res = run_cont(res,nres,cp,ic,kp,sd,fp,sp,opt)
%RUN_CONT performs the continuation process to find solutions along a branch
%for a given number of steps.
%
%INPUTS:    res         structure:  Initialised arrays for results
%           nres        structure:  Initialised arrays for results of a single 
%                                   Newton iteration
%           cp          structure:  Parameters for continuation
%           ic          structure:  Starting state for the continuation process
%           kp          structure:  Parameters for Krylov methods
%           sd          structure:  Filenames and working directory 
%           fp          structure:  Fixed parameters in the system
%           sp          structure:  Parameters for computing stability
%           opt         structure:  Paramaters for optional processes
%
%OUTPUTS    res         structure: Results along the whole branch

while nres.its <= cp.Ns
    
    % Perform Newton iteration until convergence or termination for that step.
    nres.converged = 0;
    while nres.converged == 0
        nres = find_initial_iterate(nres,fp,cp);
        nres = newton_iteration(nres,kp,fp,cp,sd,opt);
    end
    
    if nres.converged == 2
        return
    end
    
    % Compute eigenvalues if needed
    if sp.find_eval == 1
        nres = comp_eval(nres,sp,fp,ic);
    end

    % Update values and store results
    [res,nres] = update_res(res,nres,cp,fp,sp,opt);
    
    nres.its = nres.its +1;
end
end

