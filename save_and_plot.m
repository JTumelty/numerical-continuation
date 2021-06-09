function  save_and_plot(res,kp,cp,sp,ic,fp,sd,opt)
% SAVE_AND_PLOT plots several plots from the data obtained including:
%   -Bifurcation diagram in (r,N) space
%   -Cost of iterations along the branch
%   -Solution profiles at the specified positions from opt.plot_profile
%   -Rate of Newton convergence
%   -Leading 4 eigenvalues
% and saves all the computed results in a mat file.
save(strcat(sd.dir,'/',sd.dataname,'.mat'),'res','sd','opt','kp','cp','fp','sp','ic')
[res.left_sn, res.right_sn] = find_sn(res,cp,sp,fp);


% Plot bifurcation diagram
clf; hold on
plot(res.r(1:find(res.r,1,'last')),res.E(1:find(res.r,1,'last')))
xlabel('r')
ylabel('N')
saveas(gcf,strcat(sd.dir,'/bifurcation_diagram.fig'))
if cp.bt == 1 || cp.mp == 1
    clf; hold on
    plot(res.r(1:find(res.r,1,'last')),res.nu(1:find(res.r,1,'last')))
    xlabel('r')
    ylabel('\nu')
    saveas(gcf,strcat(sd.dir,'/r_nu_bifurcation_diagram.fig'))
end

% Plot cost along the branch
res.meancost = mean(res.cost);
clf; hold on
plot(abs(res.s(1:find(res.r,1,'last'))),res.cost(1:find(res.r,1,'last')))
xlabel('s')
ylabel('cost')
saveas(gcf,strcat(sd.dir,'/cost.fig'))

% Plot example solution profiles
if ~isfield(opt,'plot_profile')
    for i = opt.plot_profile
        if i < length(res.s(1:find(res.r,1,'last')))
            clf; hold on
            plot(fp.x,res.u(:,i))
            xlabel('x')
            ylabel('u')
            saveas(gcf,strcat(sd.dir,'/profile_',num2str(round(res.s(i),3)),'.fig'))
        end
    end
end

% Plot free energy profiles
clf; hold on
plot(abs(res.r(1:find(res.r,1,'last'))),res.F(1:find(res.r,1,'last')))
plot([min(res.r) max(res.r)], [0 0],'--')
xlabel('r')
ylabel('F')
saveas(gcf,strcat(sd.dir,'/free_energy.fig'))

% Plots for finding the rate of convergence.
res.p = zeros(size(res.fall));
res.p_end = zeros(length(res.fall(1,:)),1);
for i = 1:length(res.fall(1,:))
    for j = 2:find(res.fall(:,i),1,'last')
        res.p(j-1,i) = log(res.fall(j,i))/log(res.fall(j-1,i));
        res.p_end(i) = res.p(j-1,i);
    end
end
res.meanp = mean(res.p_end(2:end));
clf; hold on
plot(abs(res.s(2:end)),res.p_end(2:end),'x')
xlabel('s')
ylabel('Rate of Newton convergence')
saveas(gcf,strcat(sd.dir,'/Newton_rate.fig'))

% Plot eigenvalues
if sp.find_eval == 1
clf; hold on
for i = 1:min(sp.evalno,sp.evalno)
    s_plot = abs(res.s(1:round(sp.bps/cp.ts):end));
    plot(s_plot(1:length(res.eval(i,:))),res.eval(i,:),'.')
end
xlabel('s')
ylabel('\lambda')
saveas(gcf,strcat(sd.dir,'/eval.fig'))
end

%Save data
save(strcat(sd.dir,'/',sd.dataname,'.mat'),'res','sd','opt','kp','cp','fp','sp','ic')

end

