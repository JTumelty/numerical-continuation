function nres = check_termination(nres,cp,sd,tol)
% CHECK_TERMINATION outputs results to log file and determines if 
% convergence of the Newton step has been achieved.
%
%INPUTS: Required: nres, cp, sd
%        Optional: tol
%
%OUTPUTS: nres with updated stepsize for next step in continuation process

% Update log file with Newton iteration results. 
% If 3 inputs are given then information is from initial conditions.
if nargin == 3
    fprintf(sd.fileID, '%4.0f \t %4.10f \t %4.10f \t %4.0f \t -- \t --- f\t --- \t %4.12f \t --- \n',...
           nres.its,nres.rn,nres.ds,nres.itn,nres.f);
else 
    fprintf(sd.fileID, '%4.0f \t%4.10f \t %4.10f \t %4.0f \t %4.0f \t %4.8f\t %4.8f \t %4.12f \t %4.1f \n',...
            nres.its,nres.rn,nres.ds,nres.itn,nres.flag,tol,nres.relres,nres.f,nres.cost);   
end
 
% Check if convergence is reached and end
if nres.f < cp.toln 
       fprintf(sd.fileID, 'Converged %4.0f \n\n',nres.cost);
       nres.converged = 1;
        % Vary the step-size if convergence was achieved using too many or
        % too few iterations in total. 
        if cp.vds == 1 && isfield(cp,'mincost') && isfield(cp,'maxcost')
            if nres.cost < cp.mincost && abs(nres.ds) < cp.maxds
                nres.ds = 1.2*nres.ds;
            elseif nres.cost > cp.maxcost && abs(nres.ds) > cp.minds
                nres.ds = 0.8*nres.ds;
            end
        end   
    return        
end


% Stop if reached max number of iterations.
if (nres.itn >= cp.maxitn) 
    fprintf(sd.fileID,'Max number of Newton iterations reached \n\n');
    
    % If allowed, decrease the step-size and restart the current Newton
    % iteration from the beginning.
    %
    % Otherwise terminate the program.
    if cp.vds == 1
        if abs(nres.ds) > cp.minds
            nres.ds = nres.ds/2;
            nres.restart = 1;
            nres.cost = 0;
        else
            warning('Error. \n Maximum number of Newton iterations reached.')
            nres.converged = 2;
        end
    else
        warning('Error. \n Maximum number of Newton iterations reached.')
        nres.converged = 2;
    end
    
% Stop if there is a non-zero flag in the Krylov iteration. 
elseif nres.flag ~= 0
    fprintf(sd.fileID, 'Flag non-zero \n\n');
    
        % If allowed, decrease the step-size and restart the current Newton
        % iteration from the beginning.
        %
        % Otherwise terminate the program.
        if cp.vds == 1
            if abs(nres.ds) > cp.minds
                nres.ds = nres.ds/2;
                nres.restart = 1;
                nres.cost = 0;
            else
                warning('Error. \n Reached non-zero flag in Krylov iteration.')
                nres.converged = 2;
            end
        end
     warning('Error. \n Reached non-zero flag in Krylov iteration.')
     nres.converged = 2;
end

end

