function [res,nres] = init_res(fp, cp, ic,opt)
%INIT_RES initialises the structures res and nres to hold regular results
%and results for a single continuation step respectively.
%
%INPUTS:    fp
%           cp
%           ic
%           opt
%
%OUTPUTS:   res 
%           nres

% Initialise res
size            = ceil(cp.Ns/cp.ts);

res.u           = zeros(fp.N,size);
res.E           = zeros(1,size);
res.F           = zeros(1,size);
res.FL          = zeros(1,size);
res.r           = zeros(1,size);
res.s           = zeros(1,size);
res.nu          = zeros(1,size);
res.f           = zeros(1,size);
res.fall        = zeros(cp.maxitn+1, size);

if opt.armijo == 1
    res.lall = zeros(cp.maxitn, size);
end

res.ds_all      = zeros(1,size);

res.cost        = zeros(1,size);
res.sc          = 1;
res.bp_sc       = 1;

% Initialise nres
nres.sn         = 0;
nres.un         = ic.u;
nres.rn         = ic.r;
nres.ds         = cp.ds;
nres.nun        = ic.nu;

nres.uprev      = [ic.u zeros(fp.N,1)];
nres.rprev      = [ic.r 0];
nres.sprev      = [0 0];

nres.its        = 1;
nres.itn        = 0;
nres.restart    = 0;


if cp.mp == 1
    nres.nuprev = [ic.nu 0];
end
if cp.bt == 1
    nres.vn     = ic.v;
    nres.vprev  = [ic.v zeros(fp.N,1)];
    nres.nuprev = [ic.nu 0];
    res.v       = zeros(fp.N,size);
end

end

