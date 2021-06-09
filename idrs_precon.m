function M = idrs_precon(x,P)
%IDRS_PRECON returns the action of the chosen preconditioner (L^-1 or
%dt(I-dtL)^-1) on the vector x. 
%
%INPUTS: P      structure: Parameters rn,k,N, dt (when dt(I-dtL)^-1 is used)
%        x      vector:    State that the preconditioner acts on
%
%OUTPUTS: M     vector:    Action of preconditioner on x.

% Apply the preconditioner L^-1
if ~isfield(P,'dt')
    L = P.rn - (1 - P.k.^2).^2;
    for i = 1:length(L)
        if abs(L(i)) < 1e-6
            L(i) = 1;
        end
    end

    if length(x) == P.N
        M = ifft(fft(x)./L,'symmetric');  
    else
        M = [ifft(fft(x(1:P.N))./L,'symmetric'); x(P.N+1)];
    end
% Apply the preconditioner dt(I-dtL)^-1
else
    P = (1 - P.dt*(P.rn - (1 - P.k.^2).^2))./P.dt;
    for i = 1:length(L)
        if P(i) < 1e-6  
            P(i) = 1;
        end
    end

    if length(x) == P.N
        M = ifft(fft(x)./P,'symmetric');  
    else
        M = [ifft(fft(x(1:P.N))./P,'symmetric'); x(P.N+1)];
    end
end
end