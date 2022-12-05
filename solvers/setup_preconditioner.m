function P = setup_preconditioner(kp,k,r)
    % SETUP_PRECONDITIONER initialises the preconditioner P for the linear
    % solvers of the Swift--Hohenberg equation
    % INPUTS: kp    Krylov parameters for type of preconditioner and size
    %               dt if kp.precon = 1
    %         k     vector of wavenumbers
    %         r     Current value of the bifurcation parameter
    % OUTPUTS: P    Vector P that is used as the preconditioner in spectral
    %               space
    
    if kp.precon == 0
        % Use L^-1 as preconditioner
        P = r - (1 - k.^2).^2;
    elseif kp.precon == 1
        % Use Stokes preconditioning based on taking a time-step of size dt
        P = (1 - kp.dt*(r - (1 - k.^2).^2))./kp.dt;
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
end