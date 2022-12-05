function [f,x,flag,relres,iter] = solver_bt(nres,fp,kp,tol)
        % SOLVER_BT implemente a linear solver to solve Jx = F for
        % numerical continuation of a bifurcation point in (r,nu) parameter
        % space.
        % INPUTS:   nres      Information about the current continuation step
        %           fp        Fixed parameters needed to cast the system
        %                     into a solvable form
        %           kp        Krylov parameters used in the Krylov subspace
        %                     methods
        %           tol       Tolerance for 
        % OUTPURS:  x         Vector x that solves Jx = F
        %           flag      Indicates whether the algorithm has
        %                     successfully converged (flag = 0)
        %           relres    Relative residual
        %           iter      The number of inner and outer iterations
        %                     needed for convergence.
        if nargin == 2
            f = rhs();
        elseif nargin == 4
            if kp.kmeth == 0
                [x,flag,relres,iter] = gmres(@lhs,rhs,[],tol,kp.maxit,@precon);
            elseif kp.kmeth == 1
                [x,flag,relres,iter] = bicgstab(@lhs,rhs,tol,kp.maxit,@precon);
             else
                error('Error. \n kp.kmeth not valid. Check that it takes 0 or 1.')
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

        Fp1 = ifft((nres.rn - (1 -fp.k.^2).^2).*fft(x(1:fp.N)),'symmetric') + ...
            (2*nres.nun*nres.un - 3*nres.un.^2).*x(1:fp.N) + nres.un.*x(2*fp.N+1);
        Fp2 = (2*nres.nun*nres.vn - 6*nres.un.*nres.vn).*x(1:fp.N) + ...
            ifft((nres.rn - (1 -fp.k.^2).^2).*fft(x(fp.N+1:2*fp.N)),'symmetric') + ...
            (2*nres.nun*nres.un - 3*nres.un.^2).*x(fp.N+1:2*fp.N) +nres.vn*x(2*fp.N+1);
        Fp3 = 2*nres.vn.'*x(fp.N+1:2*fp.N);
        
        if nres.its <= 2
            Fp = [Fp1; Fp2; Fp3];
        else
            Fp1 = Fp1 + nres.un.^2*x(2*fp.N+2);
            Fp2 = Fp2 + 2*nres.un.*nres.vn.*x(2*fp.N+2);
            Fp4 = oc_lhs(nres.lambda,{nres.udot.',x(1:fp.N),1},...
                                     {nres.vdot.',x(fp.N+1:2*fp.N),1},...
                                     {nres.rdot,x(2*fp.N+1),fp.w},...
                                     {nres.nudot.',x(2*fp.N+2),fp.w});
                                 
            Fp = [Fp1; Fp2; Fp3; Fp4];
        end
    end

    function F = rhs()
        % Computes the RHS of Jx = F used in Krylov iteration method. See
        % notes for further details. 

        % Newton step only changes u for its < 2.
        % F = L(un) + N(un)

        F1 =  ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.un),'symmetric') +...
            nres.nun*nres.un.^2 - nres.un.^3;
        F2 =  ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.vn),'symmetric') +...
            2*nres.nun*nres.un.*nres.vn - 3*nres.un.^2.*nres.vn;
        F3 = norm(nres.vn)^2 - 1;
            
        if nres.its <= 2
            F = [F1; F2; F3];
        else
            F4 = oc(nres.lambda,{nres.udot.', nres.un, nres.unt,1},...
                                {nres.rdot,nres.rn,nres.rnt,fp.w},...
                                {nres.vdot.', nres.vn, nres.vnt,1},...
                                {nres.nudot,nres.nun,nres.nunt,fp.w});
            F = [F1; F2; F3; F4];
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
            M = [ifft(fft(x(1:fp.N))./P,'symmetric'); ...
                    ifft(fft(x(fp.N+1:2*fp.N))./P,'symmetric'); x(2*fp.N+1)];
        else
            % Preconditioner acts like the identity on \Delta r (N+1) equation
            % when its > 2 and pseudo-arclength continuation is used.       
             M = [ifft(fft(x(1:fp.N))./P,'symmetric'); ...
                 ifft(fft(x(fp.N+1:2*fp.N))./P,'symmetric'); ...
                 x(2*fp.N+1); x(2*fp.N+2)];
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