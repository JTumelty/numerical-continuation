function nres = check_termination(nres,cp,sd,tol)
% CHECK_TERMINATION outputs results to log file and determines if 
% convergence of the Newton step has been achieved.
%
%INPUTS: Required: nres, cp, sd
%        Optional: tol
%
%OUTPUTS: nres with updated stepsize for next step in continuation process

% Update log file with Newton iteration results. 
if nargin == 3
    output_to_log(nres,sd.fileID)
else
    output_to_log(nres,sd.fileID,tol);
end
 
% Check if convergence is achieved
if nres.f < cp.toln 
    nres = converged(nres,cp);

% Stop if reached max number of iterations.
elseif (nres.itn >= cp.maxitn) 
    nres = failed_to_converge(nres,cp,'Error. \n Maximum number of Newton iterations reached.');
    
% Stop if there is a non-zero flag in the Krylov iteration. 
elseif nres.flag ~= 0
    
    nres = failed_to_converge(nres,cp,'Error. \n Reached non-zero flag in Krylov iteration.');
end

end

%-------------------------------------------------------------------------%
function nres = converged(nres,cp)
    % CONVERGED details the processes when the solution has converged
    
    nres.converged = 1;
    % Vary the step-size if convergence was achieved using too many or
    % too few iterations in total. 
    if cp.vds == 1 && isfield(cp,'mincost') && isfield(cp,'maxcost')
        nres.ds = vary_step_size(cp,nres.ds, nres.cost);
    end
end

function ds = vary_step_size(cp, ds, cost)
    % VARY_STEP_SIZE allows the step size to vary over the continuation
    % process depending upon the cost of the Newton iterations
    
    if cost < cp.mincost && abs(ds) < cp.maxds
        % Increase the stepsize if cost is below minimum value
        ds = 1.2*ds;
    elseif cost > cp.maxcost && abs(ds) > cp.minds
        % Decrease the stepsize if cost is above the maximum
        ds = 0.8*ds;
    end
end

function nres = failed_to_converge(nres,cp,message)
    % REACHED_MAXITN details the steps to perform if the solution has not
    % converged and has reached the maximum number of iterations
    
    % If allowed, decrease the step-size and restart the current Newton
    % iteration from the beginning.
    if cp.vds == 1 &&  abs(nres.ds) > cp.minds
        nres.ds = nres.ds/2;
        nres.restart = 1;
        nres.cost = 0;
    else
        % Otherwise terminate the program.
        warning(message)
        nres.converged = 2;
    end
end
