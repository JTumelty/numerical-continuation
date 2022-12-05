function [f,x,flag,relres,iter] = solver_standard(nres,fp,kp,tol)
    % SOLVER_STANDARD implemente a linear solver to solve Jx = F for
    % numerical continuation of the standard Swift--Hohenberg equation with 
    % the parameter r as a bifurcation parameter.
    % space.
    % INPUTS:   nres      Information about the current continuation step
    %           fp        Fixed parameters needed to cast the system
    %                     into a solvable form
    %
    % Optional:
    %           kp        Krylov parameters used in the Krylov subspace
    %                     methods
    %           tol       Tolerance to use in the algorithm
    % OUTPURS:  x         Vector x that solves Jx = F
    %           flag      Indicates whether the algorithm has
    %                     successfully converged (flag = 0)
    %           relres    Relative residual
    %           iter      The number of inner and outer iterations
    %                     needed for convergence.
    
    if nargin == 2
        % With two input arguments compute the evaluation of the rhs of the system 
        f = rhs();
    elseif nargin == 4
        % With four input arguments, solve the linear system LHS x = RHS
        if kp.kmeth == 0
            [x,flag,relres,iter] = gmres(@lhs,rhs,[],tol,kp.maxit,@precon);
        elseif kp.kmeth == 1
            [x,flag,relres,iter] = bicgstab(@lhs,rhs,tol,kp.maxit,@precon);
        else
            error('Error. \n kp.kmeth not valid. Check that it takes 0, or 1.')
        end
        f = NaN;
    end

%-------------------------------------------------------------------------%
%% Solver functions for the left and right-hand sides in each Newton step 
%   when computing standard SH equation
%-------------------------------------------------------------------------%
    function Fp = lhs(x)
        % Creates a function handle for the Jacobian acting on a vector x
        % to solve Jx = F. 

        % if nres.its <= 2
        % x = du
        % J du = L du + N_u du as r is fixed for first two steps.

        Fp = ifft((nres.rn - (1 - fp.k.^2).^2).*fft(x(1:fp.N)),'symmetric')...
            + 2*nres.nun*nres.un.*x(1:fp.N) - 3*nres.un.^2.*x(1:fp.N);
        
        if nres.its > 2
            % For its > 2, 
            % J[du;dr] = [L du + N_u du + udr; lambda (udot du + w rdot dr)]
            Fp = [Fp+nres.un.*x(fp.N+1); ...
                  oc_lhs(nres.lambda,...
                                            {nres.udot.',x(1:fp.N),1},...
                                            {nres.rdot,x(fp.N+1),fp.w})];
        end
    end

    function F = rhs()
        % Computes the RHS of Jx = F used in Krylov iteration method. See
        % notes for further details. 

        % Newton step only changes u for its < 2.
        % F = L(un) + N(un)

        F = ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.un),'symmetric')...
            + nres.nun*nres.un.^2 - nres.un.^3;
        if nres.its > 2
            % Newton step changes u and r for its > 2. 
            %
            % Additional condition is an orthogonality condition imposing
            % that the difference between the Newton iterate (r,u) and initial
            % iterate (r*,u*) is orthogonal to extension from previous iterates.
            %
            % F = [L(u)+N(u); oc(r,u)]
            % where oc = lambda *(udot(u - u*) + w rdot(r - r*))
            F = [F; oc(nres.lambda,...
                                    {nres.udot.', nres.un, nres.unt,1},...
                                    {nres.rdot,nres.rn,nres.rnt,fp.w})];
        end
    end
%-------------------------------------------------------------------------%
%% Preconditioner
%-------------------------------------------------------------------------%
function M = precon(x)
        % Computes the action of either the inverse linear operator or
        % dt(I-dtL)^-1 as a preconditioner for the system Jx = F. 
        
        % Setup preconditioner for PDE
        P = setup_preconditioner(kp,fp.k,nres.rn);

        % Return the action of the preconditioner on a vector x.
        if nres.its <= 2
           M = ifft(fft(x)./P,'symmetric');  
        else
            % Preconditioner acts like the identity on \Delta r (N+1) equation
            % when its > 2 and pseudo-arclength continuation is used. 
            M = [ifft(fft(x(1:fp.N))./P,'symmetric'); x(fp.N+1)];
        end
end
end
%-------------------------------------------------------------------------%
%% Functions for solving the orthogonality condition
%-------------------------------------------------------------------------%
function F = oc(lambda,varargin)
    % lambda is a weighting factor, adds contributions like
    % udot.'(u - unt)*w  where w is a weighting.
    F = 0;
    for i = 1:nargin-1
        F = F + varargin{i}{1}*(varargin{i}{2} - varargin{i}{3})*varargin{i}{4};
    end
    F = lambda*F;
end

function F = oc_lhs(lambda,varargin)
    % Adds contributions udot.*delta u.*w to the above
    F = 0;
    for i = 1:nargin-1
        F = F + varargin{i}{1}*varargin{i}{2}*varargin{i}{3};
    end
    F = lambda*F;
end