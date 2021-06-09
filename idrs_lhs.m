function Fp = idrs_lhs(x,A)
%IDRS_LHS evaluates the action of the LHS of the linear system on x when
%IDR(s) is the chosen Krylov method (see functions lhs and rhs in
%newton_iteration for more details)
%
% A is a structure containing the field 'name' linking to this function, the
% fixed fields:       'N','k','w'
% and updated fields: 'rn','un','rdot','udot','lambda'
% for that Newton step. 

if length(x) == A.N
    Fp = ifft((A.rn - (1 - A.k.^2).^2).*fft(x),'symmetric')...
        + 2*A.nu*A.un.*x - 3*A.un.^2.*x;
else
    Fp = [(ifft((A.rn - (1 - A.k.^2).^2).*fft(x(1:A.N)),'symmetric')+...
            2*A.nu*A.un.*x(1:A.N)-3*A.un.^2.*x(1:A.N))+...
            A.un.*x(A.N+1); ...
            A.lambda*(A.udot.'*x(1:A.N)+A.w^0.5*A.rdot*x(A.N+1))];
end

end

