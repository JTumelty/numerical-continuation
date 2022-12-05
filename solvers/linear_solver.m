function  [x,flag,relres,iter] = linear_solver(nres,fp,kp,tol,cp)
    %LINEAR_SOLVER acts as an interface between the Newton iteration step 
    % and the numerical solver routines. 
    % This functions determines which system of equations is being solved 
    % and passes relevant parameters to the relevant solver.
    % INPUTS:   nres      Information about the current continuation step
    %           fp        Fixed parameters needed to cast the system
    %                     into a solvable form
    %           kp        Krylov parameters used in the Krylov subspace
    %                     methods
    %           tol       Specified tolerance for Krylov subspace methods
    %                     to reach convergence
    %           cp        Continuation parameters to determine which system
    %                     of equations need to be solved.
    % OUTPURS:  x         Vector x that solves Jx = F
    %           flag      Indicates whether the algorithm has
    %                     successfully converged (flag = 0)
    %           relres    Relative residual
    %           iter      The number of inner and outer iterations
    %                     needed for convergence.
    
    if cp.mp == 0 && cp.bt == 0
        % Standard numerical continuation
        [~,x,flag,relres,iter] = solver_standard(nres,fp,kp,tol);
    elseif cp.mp == 1 && cp.bt == 0
        % Track Maxwell point
        [~,x,flag,relres,iter] = solver_mp(nres,fp,kp,tol);
    elseif cp.bt == 1 && cp.mp == 0
        % Track a bifurcation point
        [~,x,flag,relres,iter] = solver_bt(nres,fp,kp,tol);
    else
        error('Error. Could not determine which system to solve.')
    end
    
end
