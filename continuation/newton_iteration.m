function nres = newton_iteration(nres,all_param)
% NEWTON_ITERATION starts from a given initial state for a Newton iteration
% performs a Newton-Krylov iterative process until the tolerance for Newton
% iteration is met. 
   
    % Get function handle for RHS evaluation
    if all_param.cp.mp == 0 && all_param.cp.bt == 0
        rhs = @(n, fp) solver_standard(n,fp);
    elseif all_param.cp.mp == 1 && all_param.cp.bt == 0
        rhs = @(n, fp) solver_mp(n,fp);
    elseif all_param.cp.mp == 0 && all_param.cp.bt == 1
        rhs = @(n,fp) solver_bt(n,fp);
    end

    % Update nres based on the initial guess for the solution
    nres = update_nres(nres,1,all_param.fp);
       
    % Check if initial state satisfies convergence criterion. 
    nres = check_termination(nres,all_param.cp,all_param.sd);

    % Find corrections via Newton methods until convergence or a restart
    % with a smaller step-size is needed. 
    while nres.restart == 0 && nres.converged == 0
        nres.itn = nres.itn + 1;
         
        % Set tolerance for Newton iteration
        tol = set_tol(all_param.op,nres.f,all_param.kp);

        % Solve the linear system Jx = F using chosen Krylov subspace 
        % method from gmres and bicgstab.
        [x,nres.flag,nres.relres,iter] = linear_solver(nres,all_param.fp,...
                                                all_param.kp,tol,...
                                                all_param.cp);
        nres.cost = compute_cost(nres.cost,iter,all_param.kp.kmeth);
        
        % Make the corrections at this Newton iteration
        nres = update_rn_un(nres,x,all_param.fp, all_param.cp);
        
        % Update nres with this corrected iteration
        nres = update_nres(nres,nres.itn+1,all_param.fp);
        
        % Check whether convergence has been achieved
        nres = check_termination(nres,all_param.cp,all_param.sd,tol);

    end % of Newton iteration
%-------------------------------------------------------------------------%
% END of main code
%-------------------------------------------------------------------------%

    function nres = update_nres(nres,i,fp)
        % UPDATE_NRES updates the information about the current Netwon iterate
        % by computing variational functions and the residual

        nres.f = norm(rhs(nres,fp));
        nres.fall(i) = nres.f;
        [nres.F, nres.FL] = compute_variational_functional(nres,fp);
    end
end


%-------------------------------------------------------------------------%
% Update f and make corrections to un, rn, etc
%-------------------------------------------------------------------------%
function nres = update_rn_un(nres,x,fp, cp)
    % Update current un and rn after this Newton step and find RHS.

    if cp.mp == 0 && cp.bt == 0
        update_rn_un_standard();
    elseif cp.mp == 1 && cp.bt == 0
        update_rn_un_mp();
    elseif cp.mp == 0 && cp.bt == 1
       update_rn_un_bt();
    end
    
%-------------------------------------------------------------------------% 
    function update_rn_un_standard()
        nres.un = nres.un - x(1:fp.N); 
        if nres.its > 2
            nres.rn = nres.rn - x(fp.N+1);
        end
    end

    function update_rn_un_mp()
        nres.un = nres.un - x(1:fp.N); 
        nres.rn = nres.rn - x(fp.N+1);
        if nres.its > 2
            nres.nun = nres.nun - x(fp.N+2);
        end
    end

    function update_rn_un_bt()
        nres.un = nres.un - x(1:fp.N);
        nres.vn = nres.vn - x(fp.N+1:2*fp.N);
        nres.rn = nres.rn - x(2*fp.N+1);

        if nres.its > 2
            nres.nun = nres.nun - x(2*fp.N+2);
        end
    end
end


%-------------------------------------------------------------------------%
% Compute the cost for the Krylov iteration
%-------------------------------------------------------------------------%
function cost = compute_cost(cost,iter,kmeth)
% Compute cost for each Netwon iteration, i.e., the number of times the
% linear system needs to be solved
    if kmeth == 0
        cost = cost + (iter(1)-1)*40 + iter(2);
    elseif kmeth == 1
        cost = cost + 2*iter;
    end
end
%-------------------------------------------------------------------------%
% Compute variational functional
%-------------------------------------------------------------------------%
function [F,FL] = compute_variational_functional(nres,fp)
    % Computes the free energy of the steady state
    
    F = comp_F();
    FL = comp_FL();
    
    function F = comp_F()
        L = 2*pi*fp.Lx;
        Lu = ifft((nres.rn - (1 - fp.k.^2).^2).*fft(nres.un),'symmetric');
        F = -sum(0.5*nres.un.*Lu + nres.nun*nres.un.^3/3 - nres.un.^4/4)/(fp.N*L);     
    end

    function FL = comp_FL()
        FL = -1/(2*pi^2*fp.Lx^3*fp.N)*sum(nres.un.*ifft((fp.Lx*fp.k).^2.*fft(nres.un),'symmetric') ...
            + 1/fp.Lx^2*nres.un.*ifft((fp.Lx*fp.k).^4.*fft(nres.un),'symmetric'));
    end
end
%-------------------------------------------------------------------------%
% Set tolerance 
%-------------------------------------------------------------------------%
function tol = set_tol(op,f,kp)
    % Set a combined relative-absolute error tolerance if specified,
    % otherwise use given relative tolerance. 
    if op.combtol == 1
        tol = min(0.1,op.tr + op.ta/f);
    else
        tol = kp.tolk;
    end
end
