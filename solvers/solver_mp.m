function [f,x,flag,relres,iter] = solver_mp(nres,fp,kp,tol)
    % SOLVER_MP implemente a linear solver to solve Jx = F for
    % numerical continuation of a Maxwell Point of the Swift--Hohenberg 
    % equation in (r,nu) parameter space.
    % INPUTS:   nres      Information about the current continuation step
    %           fp        Fixed parameters needed to cast the system
    %                     into a solvable form
    %           kp        Krylov parameters used in the Krylov subspace
    %                     methods
    %           tol
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
%% Solver functions for Maxwell Point
%-------------------------------------------------------------------------%
    function Fp = lhs(x)
        Fp1 = ifft((nres.rn - (1 -fp.k.^2).^2).*fft(x(1:fp.N)),'symmetric') + ...
            (2*nres.nun*nres.un - 3*nres.un.^2).*x(1:fp.N) + nres.un.*x(fp.N+1);
        
        L = 2*pi*fp.Lx;
        Lu = ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.un),'symmetric');
        Fp2 = -sum(x(1:fp.N).*(Lu + nres.nun*nres.un.^2 - nres.un.^3) + ...
            1/2*nres.un.^2*x(fp.N+1))/(fp.N*L);
        
        if nres.its <= 2
            Fp = [Fp1; Fp2];
        else
            Fp1 = Fp1 + nres.un.^2*x(fp.N+2);
            Fp2 = Fp2 - sum(nres.un.^3)/(3*fp.N*L)*x(fp.N+2);
            Fp3 = oc_lhs(nres.lambda,...
                                            {nres.udot.',x(1:fp.N),1},...
                                            {nres.rdot,x(fp.N+1),fp.w},...
                                            {nres.nudot.',x(fp.N+2),fp.w});
            
            Fp = [Fp1; Fp2; Fp3];
        end
    end

    function F = rhs()
            F1 =  ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.un),'symmetric') +...
                nres.nun*nres.un.^2 - nres.un.^3;
            F2 = comp_F();
            
        if nres.its <= 2
            F = [F1; F2];
        else
            F3 = oc(nres.lambda,...
                                            {nres.udot.', nres.un, nres.unt,1},...
                                            {nres.rdot,nres.rn,nres.rnt,fp.w},...
                                            {nres.nudot,nres.nun,nres.nunt,fp.w});
            F = [F1; F2; F3];
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
            M = [ifft(fft(x(1:fp.N))./P,'symmetric'); x(fp.N+1)];
        else
            % Preconditioner acts like the identity on \Delta r (N+1) equation
            % when its > 2 and pseudo-arclength continuation is used. 
            M = [ifft(fft(x(1:fp.N))./P,'symmetric'); x(fp.N+1); x(fp.N+2)];
        end       
    end


%-------------------------------------------------------------------------%
%% Compute free energy
%-------------------------------------------------------------------------%
    function F = comp_F()
        L = 2*pi*fp.Lx;
        Lu = ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.un),'symmetric');
        F = -sum(0.5*nres.un.*Lu + nres.nun*nres.un.^3/3 - nres.un.^4/4)/(fp.N*L);     
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