function nres = comp_eval(nres,sp,fp,ic)
%COMP_EVAL computes the leading eigenvalues and eigenvectors of the
%Jacobian, J, at specified points. 
% One of the time-stepping schemes: semi-implicit Euler, ETDRK2 or
% ETDRK4 (default) is chosen to find an approximation to the exponential of 
% the Jacobian,e^Jdt by taking a single time-step of size sp.Sdt 
% (default = 0.01 (ETDRK4),0.001 (ETDRK2), 0.0001 (Euler))
%
% eigs is used to find the (sp.evalno) leading eigenvalues of this operator
% and uses the following parameters:
% Tolerance: sp.eigstol (default 10^-8), 
% Maximum number of iterations: sp.eigsmaxit (default = 800),
% Krylov subspace dimension: sp.K (default = 20 (ETDRK4), 30 (ETDRK2), 40 (Euler))
%
% Eigenvalues of J are found from these. 

if mod(nres.its-1,sp.bps) == 0
    
    if sp.Smeth == 0
        [exp_ch,exp_ch2,f1,f2,f3,exp_ch_diff] = comp_coeff_ETDRK4();
        [V,D] = eigs(@act_etdrk4,fp.N,sp.evalno,'largestabs','Tolerance',...
            sp.eigstol,'MaxIterations',sp.eigsmaxit,'SubspaceDimension',sp.K);
    elseif sp.Smeth == 1
        [exp_ch,f,exp_ch_diff] = comp_coeff_ETDRK2();
        [V,D] = eigs(@act_etdrk2,fp.N,sp.evalno,'largestabs','Tolerance',...
            sp.eigstol,'MaxIterations',sp.eigsmaxit,'SubspaceDimension',sp.K);
    elseif sp.Smeth == 2
        [V,D] = eigs(@explhs,fp.N,sp.evalno,'largestabs','Tolerance',...
            sp.eigstol,'MaxIterations',sp.eigsmaxit,'SubspaceDimension',sp.K);
    end
    
    % If the eigenvalues of e^Jdt are mu, the eigenvalues of J are 1/dt log(mu) 
    nres.eval = 1/sp.Sdt*log(diag(D));
    nres.evec = reshape(V,[],1);
end
%-------------------------------------------------------------------------%

function F_hat = NL_F(yhat)
    %Action of N_u on a state y= ifft(yhat). 
    y       = ifft(yhat,'symmetric');
    F_hat   = fft((2*ic.nu*nres.un - 3*nres.un.^2).*y);
end

%-------------------------------------------------------------------------%
% ETDRK4
%-------------------------------------------------------------------------%
        
function [exp_ch,exp_ch2,f1,f2,f3,exp_ch_diff] = comp_coeff_ETDRK4()
    % Find the coefficients needed for ETDRK4, using a Taylor series to
    % evaluate them if dt(r - (1-k^2)^2) is small. 
    ch          = sp.Sdt*(nres.rn - (1 - fp.k.^2).^2);
    exp_ch      = exp(ch);
    exp_ch2     = exp(ch/2);

    f1           = sp.Sdt*(-4 - ch           + exp_ch.*(4 - 3*ch + ch.^2))./ch.^3;
    f2           = sp.Sdt*( 2 + ch           + exp_ch.*(-2 + ch))./ch.^3;
    f3           = sp.Sdt*(-4 - 3*ch - ch.^2 + exp_ch.*(4 - ch))./ch.^3;
    exp_ch_diff  = sp.Sdt*(exp_ch2-1)./ch;

    % Replace coefficients ch below threshold with Taylor expansion
    for chi = 1:length(ch)
        if abs(ch(chi)) < 1e-2 || ch(chi) == Inf
            exp_ch_diff(chi) = sp.Sdt*(1/2 + 1/8*ch(chi)  + 1/48*ch(chi).^2  + 1/384*ch(chi).^3 + 1/3840*ch(chi).^4);
            f1(chi)          = sp.Sdt*(1/6 + 1/6*ch(chi)  + 3/40*ch(chi).^2  + 1/45*ch(chi).^3);
            f2(chi)          = sp.Sdt*(1/6 + 1/12*ch(chi) + 1/40*ch(chi).^2  + 1/180*ch(chi).^3);
            f3(chi)          = sp.Sdt*(1/6             - 1/120*ch(chi).^2  - 1/360*ch(chi).^3);
        end
    end
end

function v = act_etdrk4(x)
    % Apply an ETDRK4 scheme to time-step the state x which
    % evolves under the action of the Jacobian. 
        x_hat   = fft(x);
        Fu      = NL_F(x_hat);
        
        a_hat   = x_hat.*exp_ch2 + Fu.*exp_ch_diff;
        Fa      = NL_F( a_hat);
        
        b_hat   = x_hat.*exp_ch2 + Fa.*exp_ch_diff;
        Fb      = NL_F(b_hat);
        
        c_hat   = a_hat.*exp_ch2 + (2*Fb - Fu).*exp_ch_diff;
        Fc      = NL_F(c_hat);
        
        u_hat   = x_hat.*exp_ch + (Fu.*f1 + 2*(Fa + Fb).*f2 + Fc.*f3);
                                
        v       = ifft(u_hat,'symmetric');
end
%-------------------------------------------------------------------------%
% ETDRK2
%-------------------------------------------------------------------------%
function [exp_ch,f,exp_ch_diff] = comp_coeff_ETDRK2()
    % Find the coefficients needed for ETDRK2, using a Taylor series to
    % evaluate them if dt(r - (1-k^2)^2) is small. 
    ch = sp.Sdt*(nres.rn - (1 - fp.k.^2).^2);  
    exp_ch = exp(ch);
    exp_ch_diff = sp.Sdt*(exp_ch-1)./ch;
    f = sp.Sdt*(exp_ch - 1 - ch)./ch.^2;

    % Replace coefficients ch below threshold with Taylor expansion
    for i = 1:length(ch)      
        if abs(ch(i))< 1e-2 || isinf(ch(i))
            exp_ch_diff(i) = sp.Sdt*(1 + 1/2*ch(i) + 1/6*ch(i)^2 + 1/24*ch(i)^3 + 1/120*ch(i)^4 + 1/720*ch(i)^5 +1/4940*ch(i)^6);
            f(i) = sp.Sdt*(1/2 + 1/6*ch(i) + 1/24*ch(i)^2 + 1/120*ch(i)^3+1/720*ch(i)^4 + 1/4940*ch(i)^5+ 1/(4940*8)*ch(i)^6);
        end
    end
end

function v = act_etdrk2(x)
    % Apply a ETDRK2 scheme to time-step the state x which
    % evolves under the action of the Jacobian. 
    x_hat   = fft(x);
    F       = NL_F(x_hat);
    a_hat   = x_hat.*exp_ch + F.*exp_ch_diff;
    v       = ifft(a_hat + f.*(NL_F(a_hat) - F),'symmetric');
end

%-------------------------------------------------------------------------%
% Semi-implicit Euler scheme
%-------------------------------------------------------------------------%

function Fp = explhs(x)
    % Apply a semi-implicit Euler scheme to time-step the state x which
    % evolves under the action of the Jacobian. 
        Fp = ifft(fft(x + sp.Sdt*(2*ic.nu*nres.un.*x -3*nres.un.^2.*x))./...
            (1 - sp.Sdt*(nres.rn - (1 -fp.k.^2).^2)),'symmetric');
end

end


