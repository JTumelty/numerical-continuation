function  all_res = run_cont(all_res,all_param)
%RUN_CONT performs the continuation process to find solutions along a branch
%for a given number of steps.
%
%INPUTS:    all_res     structure:  Initialised arrays for results
%           all_param   structure:  Parameters for continuation
%
%OUTPUTS    all_res         structure: Results along the whole branch


% Continue while below maximum number of continuation steps 
while all_res.nres.its <= all_param.cp.Ns
    
    % Perform Newton iteration until convergence or termination for that step.
    all_res.nres.converged = 0;
    while all_res.nres.converged == 0
        
        % Display status of iteration on screen
        disp(strcat('Iteration: ',num2str(all_res.nres.its), ', r=',...
            num2str(all_res.nres.rn, '%6.6f')));
        
        % Find an initial guess for prediction stage
        all_res.nres = find_initial_iterate(all_res.nres,...
                                            all_param.fp,...
                                            all_param.cp);
        
        % Perform Newton iterations on the initial guess to achieve
        % convergence
        all_res.nres = newton_iteration(all_res.nres,all_param);
    end
    
    % End program if convergence is not achieved
    if all_res.nres.converged == 2
        save(strcat(all_param.sd.dir,'/',all_param.sd.dataname),'all_res');
        error('Error. Convergence not achieved. Exiting program.')
    end
    
    % Compute eigenvalues if needed
    if all_param.sp.find_eval == 1
        all_res.nres = comp_eval(all_res.nres,all_param.sp,all_param.fp,all_param.ic);
    end
        
    % Update values and store results
    all_res = update_all_res(all_res,all_param);
    
    % Increase iteration count
    all_res.nres.its = all_res.nres.its +1;
end

end

