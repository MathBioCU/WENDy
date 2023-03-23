%% generate data

odes={'wendydata_Logistic_Growth.mat',...
    'wendydata_Lotka_Volterra.mat',...
    'wendydata_FitzHugh-Nagumo.mat',...
    'wendydata_Hindmarsh-Rose.mat',...
    'wendydata_biochemM1.mat'};

ode_num = 1;                       % select ODE from list above
load(odes{ode_num},'t','x','features','params','x0','true_vec','rhs_p');

subsamp = 8;                       % subsample data in time
tobs = t(1:subsamp:end); xsub = x(1:subsamp:end,:);
[M,nstates] = size(xsub);

%% noise model

noise_dist = 0;
noise_alg = 0;

%% set WENDy params

%%% set-and-forget params
wendy_snf_params;

%%% set weak integration
phifun = phifuns{1};                % defined in wendy_snf_params.m
meth = 'direct';                     % 'mtmin','FFT','direct','timefrac'
K_max = 5000;
K_min = length(true_vec);
mt_max = max((M-K_min)/2,1);
mt_min = 2;%rad_select(tobs,xobs,phifun,1,submt,0,1,2,mt_max,[]);

%%

nrs = 10.^(-6:-1);
runs = 200;
mts = 2:floor((length(tobs)-1)/2);
mtmins = zeros(length(nrs),runs);
for nn=1:length(nrs)
    noise_ratio = nrs(nn);
    errs_temp = zeros(length(mts),runs);
    for rr=1:runs
        rng(rr);
        rng_seed = rng().Seed; rng(rng_seed);
        [xobs,noise,~,sigma] = gen_noise(xsub,noise_ratio,noise_dist,noise_alg);
        mtmins(nn,rr) = rad_select(tobs,xobs,phifun,1,submt,0,1,2,mt_max,[]);
    end
end


%% run WENDy

nrs = 10.^(-6:-1);
runs = 200;
mts = 2:floor((length(tobs)-1)/2);
out = zeros(length(mts),length(nrs),runs);
mtmins = zeros(length(nrs),runs);
for nn=1:length(nrs)
    noise_ratio = nrs(nn);
    errs_temp = zeros(length(mts),runs);
    for rr=1:runs
        rng(rr);
        rng_seed = rng().Seed; rng(rng_seed);
        [xobs,noise,~,sigma] = gen_noise(xsub,noise_ratio,noise_dist,noise_alg);
        mtmins(nn,rr) = rad_select(tobs,xobs,phifun,1,submt,0,1,2,mt_max,[]);
        parfor mm=1:length(mts)
            mt_params = mts(mm);
            mt_cell = cellfun(@(x,y) [x,{y}], repmat({{phifun,meth}},length(mt_params),1),num2cell(mt_params(:)),'uni',0);
            tic;
            [~,~,~,~,~,~,errs,~,~,...
                Theta_cell,mt,x_jac,Jac_mat,G_0,b_0,RT,stdW,mseW,CovW] = ...
                wendy_fcn(...
                xobs,tobs,features,true_vec,toggle_smooth,...
                mt_cell,mt_min,mt_max,K_min,K_max,center_scheme,toggle_VVp_svd,...
                w0,optim_params,err_norm,iter_diff_tol,max_iter,diag_reg,pvalmin,check_pval_it);
            disp(toc)
            disp(errs(end))
            errs_temp(mm,rr) = errs(end);
            disp([noise_ratio rr mt_params])
        end
    end
    out(:,nn,:) = errs_temp;
end

%%

clf
avg_out = mean(out,3);
avg_mtmins = mean(mtmins,2);
fl_avg_mtmins = floor(avg_mtmins);
alpha = avg_mtmins - fl_avg_mtmins; 
y_mtmin = alpha*0;
for nn=1:length(nrs)
    ind = find(mts==fl_avg_mtmins(nn));
    y_mtmin(nn) = avg_out(ind,nn).*(1-alpha(nn))+avg_out(ind+1,nn).*alpha(nn);
end
loglog(mts,fliplr(avg_out), 'o-','linewidth',2); hold on
loglog(avg_mtmins,y_mtmin,'x','linewidth',4, 'markersize',14)
xlim([mts(1) mts(end)])
xlabel('m_t')
ylabel('E_2')
grid on
LEGS = {'$\underline{m}_t$','$\sigma_{NR}=10^{-6}$','$\sigma_{NR}=10^{-5}$','$\sigma_{NR}=10^{-4}$','$\sigma_{NR}=10^{-3}$','$\sigma_{NR}=10^{-2}$','$\sigma_{NR}=10^{-1}$'};
legend(fliplr(LEGS),'location','southwest','interpreter','latex','box','off', 'fontsize',13);
% columnlegend(2,LEGS,'location','northeast','interpreter','latex','box','off', 'fontsize',15);
set(gca,'fontsize',12,'YTick',10.^(-6:0),'XTick',[2 3 4 5 6 7 8 10 13 16 20 25 30])
