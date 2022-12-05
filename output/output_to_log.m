function output_to_log(nres,fileID,tol)
    % Outputs the information about each Newton iterate to the log file. 
    if nargin == 2
        % Foe each continuation step, the first line indicates how close
        % the inital condition is to the steady state.
        fprintf(fileID, '%4.0f \t %4.10f \t %4.10f \t %4.0f \t -- \t --- f\t --- \t %4.12f \t --- \n',...
               nres.its,nres.rn,nres.ds,nres.itn,nres.f);
    else 
        % Subsequent iterates will also output informaion about convergence
        fprintf(fileID, '%4.0f \t%4.10f \t %4.10f \t %4.0f \t %4.0f \t %4.8f\t %4.8f \t %4.12f \t %4.1f \n',...
                nres.its,nres.rn,nres.ds,nres.itn,nres.flag,tol,nres.relres,nres.f,nres.cost);   
    end
end