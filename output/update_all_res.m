function all_res = update_all_res(all_res,all_param)
% UPDATE_ALL_RES updates the structure res to include the computed results from
% the Newton iteration at the given intervals, including stability when
% computed.

% Update the field res.rE at every continuation step
all_res.res.rE(:,all_res.res.rE_sc) = [all_res.nres.rn, compute_energy(all_res.nres.un)];
all_res.res.rE_sc = all_res.res.rE_sc + 1;

% Update all results based on the current iterate at regular intervals
if mod(all_res.nres.its-1,all_param.cp.ts) == 0
    all_res.res = update_res(   all_res.res, ...
                                all_res.nres);
end

% At regular intervals cp.bps update stability results when they are computed
if all_param.sp.find_eval == 1 && mod(all_res.nres.its-1,sp.bps) == 0
    all_res.stab_res = update_stab_res( all_res.stab_res,...
                                        all_res.nres);
end

% Update nres
all_res.nres = update_nres( all_res.nres);

% Save results
if all_res.res.sc-1 == all_res.res.maxsc
    save_results(all_res,all_param);
    % Reset res and stab_res structures
    all_res = init_all_res(all_param, all_res);
end

end
%-------------------------------------------------------------------------%
%% Compute energy
%-------------------------------------------------------------------------%
function E = compute_energy(u)
    E = (sum(abs(u).^2)/size(u,1))^0.5;
end
%-------------------------------------------------------------------------%
%% Update res
%-------------------------------------------------------------------------%
function res = update_res(res, nres)
    %UPDATE_RES updates the streuctre res with results from a single
    %iterate stored in nres.
    res.s(res.sc)                   = nres.sn;
    res.ds_all(res.sc)              = nres.ds;
    res.u(:,res.sc)                 = nres.un; 
    res.E(res.sc)                   = compute_energy(nres.un);
    res.F(res.sc)                   = nres.F;
    res.FL(res.sc)                  = nres.FL;
    res.r(res.sc)                   = nres.rn;
    
    res.f(res.sc)                   = nres.f;
    res.fall(:,res.sc)              = nres.fall;
    res.cost(res.sc)                = nres.cost;

    if isfield(res,'nu') && isfield(nres,'nun')
        res.nu(res.sc)              = nres.nun;
    end
    
    if isfield(res,'v') && isfield(nres,'vn')
        res.v(:,res.sc)             = nres.vn; 
    end
    
    res.sc = res.sc + 1;
end
%-------------------------------------------------------------------------%
%% Update current iterate
%-------------------------------------------------------------------------%
function nres = update_nres(nres)
    %UPDATE_NRES updates relevant fields in the structure nres
    
    % Update the values of the two previous iterates
    nres.rprev = [nres.rprev(2) nres.rn];
    nres.sprev = [nres.sprev(2) nres.sn];
    nres.uprev = [nres.uprev(:,2) nres.un];

    if isfield(nres,'nuprev') && isfield(nres,'nun')
        nres.nuprev = [nres.nuprev(2) nres.nun];
    end
    if isfield(nres,'vprev') && isfield(nres,'vn')
        nres.vprev = [nres.vprev(:,2) nres.vn];
    end
    
    % Increment the position along the branch.
    nres.sn    = nres.sn + nres.ds;
    
    % Reset information about Newton iterations
    nres.fall = zeros(size(nres.fall));

    nres.cost = 0;
    nres.flag = 0;
    nres.relres = 0;
    nres.itn = 0;
end
%-------------------------------------------------------------------------%
%% Update stwability
%-------------------------------------------------------------------------%
function stab_res = update_stab_res(stab_res, nres)
    %UPDATE_STAB_RES Updates the stability results.
    
    stab_res.evalrE(:,stab_res.bp_sc) = [nres.rn; compute_energy(nres.un)];
    stab_res.eval(:,stab_res.bp_sc) = nres.eval;
    stab_res.evec(:,stab_res.bp_sc) = nres.evec;
    stab_res.bp_sc = stab_res.bp_sc + 1;
    
end
%-------------------------------------------------------------------------%
%% Save results
%-------------------------------------------------------------------------%
function save_results(all_res, all_param)
    %SAVE_RESULTS saves the structure all_res in a file with path:
    % sd.dir/sd.dataname_n.mat, where n is the current results number
    save(strcat(all_param.sd.dir,'/',all_param.sd.dataname,'_',...
        num2str(all_res.res.res_no),'.mat'),'all_res');
    
end