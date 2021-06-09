function [res,nres] = update_res(res,nres,cp,fp,sp,opt)
% UPDATE_RES updates the structure res to include the computed results from
% the Newton iteration at the given intervals, including stability when
% computed.

% Each iteration update previous values of r,s and u
nres.rprev = [nres.rprev(2) nres.rn];
nres.sprev = [nres.sprev(2) nres.sn];
nres.uprev = [nres.uprev(:,2) nres.un];

if cp.mp == 1
    nres.nuprev = [nres.nuprev(2) nres.nun];
elseif cp.bt == 1
    nres.nuprev = [nres.nuprev(2) nres.nun];
    nres.vprev = [nres.vprev(:,2) nres.vn];
end

% At regular intervals cp.ts update general results fields
if mod(nres.its-1,cp.ts) == 0
    res.s(res.sc)                   = nres.sn;
    res.ds_all(res.sc)              = nres.ds;
    res.u(:,res.sc)                 = nres.un; 
    res.E(res.sc)                   = (sum(abs(nres.un).^2)/fp.N)^0.5;
    res.F(res.sc)                   = nres.F;
    res.FL(res.sc)                  = nres.FL;
    res.r(res.sc)                   = nres.rn;
    res.nu(res.sc)                  = nres.nun;
    res.f(res.sc)                   = nres.f;
    res.fall(:,res.sc)              = nres.fall;
    res.cost(res.sc)                = nres.cost;
    
    if opt.armijo == 1
        res.lall(:,res.sc)    = nres.lall;
    end
    
    if cp.bt == 1
        res.v(:,res.sc)       = nres.vn; 
    end
    
    res.sc                    = res.sc + 1;
end

% At regular intervals cp.bps update stability results when they are computed
if sp.find_eval == 1
    if mod(nres.its-1,sp.bps) == 0
        res.eval(:,res.bp_sc) = nres.eval;
        res.evec(:,res.bp_sc) = nres.evec;
        res.bp_sc = res.bp_sc + 1;
    end
end

% Increment the position along the branch.
nres.sn    = nres.sn + nres.ds;
end

