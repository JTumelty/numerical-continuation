function nres = newton_iteration(nres,kp,fp,cp,sd,opt)
% NEWTON_ITERATION starts from a given initial state for a Newton iteration
% performs a Newton-Krylov iterative process until the tolerance for Newton
% iteration is met. 

    % Reset nres results for initial state
    nres.fall = zeros(cp.maxitn+1,1);
    if opt.armijo == 1
        nres.lall = zeros(cp.maxitn,1);
    end
    if cp.bt == 0
        nres.f = norm(rhs);
    else
        nres.f = norm(rhs_bt);
    end
    nres.fall(1) = nres.f;
    nres.cost = 0;
    nres.flag = 0;
    nres.relres = 0;
    nres.F = comp_F();
    nres.FL = comp_FL();
    
    % Check if initial state satisfies convergence criterion. 
    nres = check_termination(nres,cp,sd);

    % Find corrections via Newton methods until convergence or a restart
    % with a smaller step-size is needed. 
    while nres.restart == 0 && nres.converged == 0
        nres.itn = nres.itn + 1;
         
        % Set a combined relative-absolute error tolerance if specified,
        % otherwise use given relative tolerance. 
        if opt.combtol == 1
            tol = set_tol();
        else
            tol = kp.tolk;
        end

        % Solve the linear system Jx = F using chosen Krylov subspace 
        % method from gmres, bicgstab and IDR(s).
        
        % Tracking the Maxwell Point
        if cp.mp == 1
            if kp.kmeth == 0
                [x,nres.flag,nres.relres,iter] = gmres(@lhs_mp,rhs_mp,[],tol,kp.maxit,@precon);
            elseif kp.kmeth == 1
                [x,nres.flag,nres.relres,iter] = bicgstab(@lhs_mp,rhs_mp,tol,kp.maxit,@precon);
            elseif kp.kmeth == 2
                % TODO
                %update_idrs_struct();
                %[x,nres.flag,nres.relres,iter] = idrs(kp.iLHS,rhs_bt,kp.s,tol,kp.maxit,kp.iP);
            else
                error('Error. \n kp.kmeth not valid. Check that it takes 0, 1 or 2.')
            end   
        % Standard bifurcation diagram
        elseif cp.bt == 0
            
            if kp.kmeth == 0
                [x,nres.flag,nres.relres,iter] = gmres(@lhs,rhs,[],tol,kp.maxit,@precon);
            elseif kp.kmeth == 1
                [x,nres.flag,nres.relres,iter] = bicgstab(@lhs,rhs,tol,kp.maxit,@precon);
            elseif kp.kmeth == 2
                update_idrs_struct();
                [x,nres.flag,nres.relres,iter] = idrs(kp.iLHS,rhs,kp.s,tol,kp.maxit,kp.iP);
            else
                error('Error. \n kp.kmeth not valid. Check that it takes 0, 1 or 2.')
            end
            
        % (r, nu) bifurcation tracking
        elseif cp.bt == 1
            if kp.kmeth == 0
                [x,nres.flag,nres.relres,iter] = gmres(@lhs_bt,rhs_bt,[],tol,kp.maxit,@precon);
            elseif kp.kmeth == 1
                [x,nres.flag,nres.relres,iter] = bicgstab(@lhs_bt,rhs_bt,tol,kp.maxit,@precon);
            elseif kp.kmeth == 2
                % TODO
                %update_idrs_struct();
                %[x,nres.flag,nres.relres,iter] = idrs(kp.iLHS,rhs_bt,kp.s,tol,kp.maxit,kp.iP);
            else
                error('Error. \n kp.kmeth not valid. Check that it takes 0, 1 or 2.')
            end
        end
        
        % Make the corrections at this Newton iteration
        update_rn_un()
        comp_cost(iter)     
        nres.F = comp_F();
        nres.FL = comp_FL();
        nres = check_termination(nres,cp,sd,tol);

    end % of Newton iteration
%-------------------------------------------------------------------------%
% END of main code
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Set tolerance 
%-------------------------------------------------------------------------%
    function tol = set_tol()
        tol = min(0.1,opt.tr + opt.ta/nres.f);
    end

%-------------------------------------------------------------------------%
% Preconditioner
%-------------------------------------------------------------------------%
    function M = precon(x)
        % Computes the action of either the inverse linear operator or
        % dt(I-dtL)^-1 as a preconditioner for the system Jx = F. 
        
        if kp.precon == 0
            % Use L^-1 as preconditioner
            P = nres.rn - (1 - fp.k.^2).^2;
            
        elseif kp.precon == 1
            P = (1 - kp.dt*(nres.rn - (1 - fp.k.^2).^2))./kp.dt;

        else 
            error('Error. \n kp.precon not valid. Check that it takes 0 or 1.')
        end
        % The limit 1e-6 in the preconditioners is chosen to prevent the
        % preconditioner being (nearly) singular. (Optimal value not tested)        
        for i = 1:length(P)
                if abs(P(i)) < 1e-6
                    P(i) = 1;
                end
        end

        % Return the action of the preconditioner on a vector x.
        if nres.its <= 2
            if cp.bt == 0
                if cp.mp == 1
                    M = [ifft(fft(x(1:fp.N))./P,'symmetric'); x(fp.N+1)];
                else
                    M = ifft(fft(x)./P,'symmetric');  
                end
            elseif cp.bt == 1
                M = [ifft(fft(x(1:fp.N))./P,'symmetric'); ...
                    ifft(fft(x(fp.N+1:2*fp.N))./P,'symmetric'); x(2*fp.N+1)];
            end
        else
            % Preconditioner acts like the identity on \Delta r (N+1) equation
            % when its > 2 and pseudo-arclength continuation is used. 
            if cp.bt == 0
                if cp.mp == 1
                    M = [ifft(fft(x(1:fp.N))./P,'symmetric'); x(fp.N+1); x(fp.N+2)];
                else
                    M = [ifft(fft(x(1:fp.N))./P,'symmetric'); x(fp.N+1)];
                end
            elseif cp.bt == 1
                 M = [ifft(fft(x(1:fp.N))./P,'symmetric'); ...
                     ifft(fft(x(fp.N+1:2*fp.N))./P,'symmetric'); ...
                     x(2*fp.N+1); x(2*fp.N+2)];
            end
        end
        
    end
%-------------------------------------------------------------------------%
%  Functions for the left and right-hand sides in each Newton step
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
                  oc_lhs(nres.lambda,{nres.udot.',x(1:fp.N),1},...
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
            F = [F; oc(nres.lambda,{nres.udot.', nres.un, nres.unt,1},...
                                   {nres.rdot,nres.rn,nres.rnt,fp.w})];
        end
    end
%-------------------------------------------------------------------------%
% Maxwell point
%-------------------------------------------------------------------------%
    function Fp = lhs_mp(x)
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
            Fp3 = oc_lhs(nres.lambda,{nres.udot.',x(1:fp.N),1},...
                                     {nres.rdot,x(fp.N+1),fp.w},...
                                     {nres.nudot.',x(fp.N+2),fp.w});
            
            Fp = [Fp1; Fp2; Fp3];
        end
    end

    function F = rhs_mp()
            F1 =  ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.un),'symmetric') +...
                nres.nun*nres.un.^2 - nres.un.^3;
            F2 = comp_F();
            
        if nres.its <= 2
            F = [F1; F2];
        else
            F3 = oc(nres.lambda,{nres.udot.', nres.un, nres.unt,1},...
                                {nres.rdot,nres.rn,nres.rnt,fp.w},...
                                {nres.nudot,nres.nun,nres.nunt,fp.w});
            F = [F1; F2; F3];
        end
        
    end
%-------------------------------------------------------------------------%
% Bifurcation tracking
%-------------------------------------------------------------------------%
    function Fp = lhs_bt(x)
        
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

    function F = rhs_bt()
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
% Orthogonality condition
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
        % Adds contributions udot.*delta u.*w
        F = 0;
        for i = 1:nargin-1
            F = F + varargin{i}{1}*varargin{i}{2}*varargin{i}{3};
        end
        F = lambda*F;
    end
%-------------------------------------------------------------------------%
% Armijo's Rule
%-------------------------------------------------------------------------%
    function armijo(x,f_ind)
        % Use of Armijo's rule may help convergence of Newton iteration by
        % decreasing the step-size taken so that F always decreases. 
        
        % Keep a record of starting un and rn.
        nres.l = 0;
        ustart = nres.un;
        rstart = nres.rn;
        
        % Make initial correction
        nres.un = nres.un - x(1:fp.N);
        if length(x) == fp.N+1
            nres.rn = nres.rn - x(fp.N+1);
        end
        nres.f = norm(rhs());  
        
        % Decrease step-size until convergence is achieved or maximum
        % number of steps (opt.maxsteps) is reached. 
        while nres.f > (1 - opt.aalpha*2^nres.l)*f_ind && nres.l > opt.maxsteps
            nres.l = nres.l-1;
            nres.un = ustart - 2^nres.l*x(1:fp.N);
            if length(x) == fp.N+1
                nres.rn = rstart - 2^nres.l*x(fp.N+1);
            end
            nres.f = norm(rhs());
        end
        
        % Record number of decreases in step size that have been taken. 
        nres.lall(nres.itn) = nres.l+1;
    end

%-------------------------------------------------------------------------%
% Functions if IDR(s) is used
%-------------------------------------------------------------------------%
    function  update_idrs_struct()
        % Update IDR(s) structures for un and rn before Krylov-iteration. 
        kp.iLHS.un = nres.un;
        kp.iLHS.rn = nres.rn;
        if nres.its > 2
            kp.iLHS.rdot = nres.rdot;
            kp.iLHS.lambda = nres.lambda;
            kp.iLHS.udot = nres.udot;
        end
        kp.iP.rn = nres.rn;
    end
%-------------------------------------------------------------------------%
% Update f and make corrections to un, rn, etc
%-------------------------------------------------------------------------%
    function update_rn_un()
        % Update current un and rn after this Newton step and find RHS.
        %            
        % Implement Armijo rule if needed.
        if cp.mp == 1
            nres.un = nres.un - x(1:fp.N); 
            nres.rn = nres.rn - x(fp.N+1);
            if nres.its > 2
                nres.nun = nres.nun - x(fp.N+2);
            end
            nres.f = norm(rhs_mp);
        elseif cp.bt == 0
            if opt.armijo == 0
                nres.un = nres.un - x(1:fp.N); 
                
                if nres.its > 2
                    nres.rn = nres.rn - x(fp.N+1);
                end
                nres.f = norm(rhs);

            elseif opt.armijo == 1
                armijo(x,nres.f); 
            else
                error('Error. \n opt.armijo not valid. Check that it takes 0 or 1.')
            end
            
        elseif cp.bt == 1
            if opt.armijo ~= 0
                error('Cannot apply armijo rule with bifurcation tracking.')
            end
            
            nres.un = nres.un - x(1:fp.N);
            nres.vn = nres.vn - x(fp.N+1:2*fp.N);
            nres.rn = nres.rn - x(2*fp.N+1);

            if nres.its > 2
                nres.nun = nres.nun - x(2*fp.N+2);
            end
           nres.f = norm(rhs_bt);
        end
        nres.fall(nres.itn+1) = nres.f;
    end
%-------------------------------------------------------------------------%
% Compute the cost for the Krylov iteration
%-------------------------------------------------------------------------%
    function comp_cost(iter)
    % Compute cost for iteration.      
        if kp.kmeth == 0
            nres.cost = nres.cost + (iter(1)-1)*40 + iter(2);
        elseif kp.kmeth == 1
            nres.cost = nres.cost + 2*iter;
        elseif kp.kmeth == 2
            nres.cost = nres.cost + iter;
        end
    end
%-------------------------------------------------------------------------%
% Compute variational functional
%-------------------------------------------------------------------------%
    function F = comp_F()
        L = 2*pi*fp.Lx;
        Lu = ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.un),'symmetric');
        F = -sum(0.5*nres.un.*Lu + nres.nun*nres.un.^3/3 - nres.un.^4/4)/(fp.N*L);     
    end

    function FL = comp_FL()
        FL = -1/(2*pi^2*fp.Lx^3*fp.N)*sum(nres.un.*ifft((fp.Lx*fp.k).^2.*fft(nres.un),'symmetric') ...
            + 1/fp.Lx^2*nres.un.*ifft((fp.Lx*fp.k).^4.*fft(nres.un),'symmetric'));
    end
%-------------------------------------------------------------------------%
end
