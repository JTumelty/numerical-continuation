function [left_sn,right_sn] = find_sn(res,cp,sp,fp)
%FIND_SN finds saddle-nodes of the branch of solutions
if sp.find_eval == 0
    left_sn = zeros(4,50);
    right_sn = zeros(4,50);
else
    left_sn = zeros(4 + fp.N,50);
    right_sn = zeros(4+fp.N,50);
end
lc = 1;
rc = 1;
for i = 2:min(length(res.r))-1
    if res.r(i) > res.r(i+1) && res.r(i) > res.r(i-1)
        right_sn(1:4,rc) = [i; res.r(i); res.E(i); res.s(i)];
        if sp.find_eval == 1
            bpind = round((i-1)/(sp.bps/cp.ts) + 1);
            if bpind < res.bp_sc
                right_sn(5:fp.N+4,rc) = find_evec(i,bpind);
            else
                break
            end
            
        end
        rc = rc + 1;
    elseif res.r(i) < res.r(i+1) && res.r(i) < res.r(i-1)
        left_sn(1:4,lc) = [i; res.r(i); res.E(i); res.s(i)];
        if sp.find_eval == 1
            bpind = round((i-1)/(sp.bps/cp.ts) + 1);
            if bpind < res.bp_sc
                left_sn(5:fp.N+4,rc) = find_evec(i,bpind);
            else
                break
            end
        end
        lc = lc + 1;
    end
end

left_sn = left_sn(:,1:lc-1);
right_sn = right_sn(:,1:rc-1);

    function v = find_evec(i,bpind)
        for j = 1:sp.evalno
            if res.evec((j-1)*fp.N+(fp.N/4+3),bpind)*res.evec(j*fp.N-(fp.N/4+3),bpind) > 0
                v = res.evec((j-1)*fp.N+1:j*fp.N,bpind);
                break
            end
        end
    end
end

