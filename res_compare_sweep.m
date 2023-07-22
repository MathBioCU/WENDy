%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

runs = 500;
numbins = 80;
xbins = linspace(-5,5,numbins+1);
xplot = (xbins(1:end-1)+xbins(2:end))/2;
pdf_true = 1/sqrt(2*pi)*exp(-xplot.^2/2);

res_wendy = zeros(runs,numbins);
res_ols = zeros(runs,numbins);

    %%% generate data
odes={'wendydata_Logistic_Growth.mat',...
    'wendydata_Lotka_Volterra.mat',...
    'wendydata_FitzHugh-Nagumo.mat',...
    'wendydata_Hindmarsh-Rose.mat',...
    'wendydata_biochemM1.mat'};

for ode_num = 1:3                       % select ODE from list above
load(odes{ode_num},'t','x','features','params','x0','true_vec','rhs_p');
% custom_ode;

%%% subsample timepoints

if ode_num==1
    subsamp = 2;                       % subsample data in time
else
    subsamp = 4;
end
tobs = t(1:subsamp:end); xsub = x(1:subsamp:end,:);
[M,nstates] = size(xsub);

%%% add noise

noise_ratio = 0.2;
noise_dist = 0;
noise_alg = 0;
wendy_snf_params;

parfor rr=1:runs

    rng(rr);
    rng_seed = rng().Seed; rng(rng_seed);
    
    [xobs,noise,~,sigma] = gen_noise(xsub,noise_ratio,noise_dist,noise_alg);
    
    %%% set weak integration
    phifun = phifuns{1};                % defined in wendy_snf_params.m
    meth = 'mtmin';                    % 'mtmin','FFT','direct','timefrac'
    mt_params = 2.^(0:3);               % see get_rad.m
    K_max = 5000;
    K_min = length(true_vec);
    mt_max = max(floor((M-1)/2)-K_min,1);
    mt_min = rad_select(tobs,xobs,phifun,1,submt,0,1,2,mt_max,[]);
    mt_cell = cellfun(@(x,y) [x,{y}], repmat({{phifun,meth}},length(mt_params),1),num2cell(mt_params(:)),'uni',0);
    
    use_true = 1;
    tic;
    [w_hat,res,res_true,res_0,res_0_true,w_hat_its,errs,V_cell,Vp_cell,...
        Theta_cell,mt,xobs,Jac_mat,G_0,b_0,RT,stdW,mseW,CovW] = wendy_fcn(xobs,tobs,features,true_vec,toggle_smooth,...
        mt_cell,mt_min,mt_max,K_min,K_max,center_scheme,toggle_VVp_svd,...
        w0,optim_params,err_norm,iter_diff_tol,max_iter,diag_reg,pvalmin,check_pval_it,use_true);
    total_time = toc;

    r_ols = b_0 - G_0*true_vec;
    r_wendy = RT \ r_ols;
        
    res_wendy(rr,:) = histcounts(r_wendy/std(r_wendy),xbins,'normalization','pdf');
    res_ols(rr,:) = histcounts(r_ols/std(r_ols),xbins,'normalization','pdf');
    disp(rr)
    
end

plot(xplot,mean(res_wendy),'r',xplot,mean(res_ols),'b',xplot,pdf_true,'g--','linewidth',2)
legend({'WENDy','OLS','N(0,1)'})
xlabel('r','fontsize',16)
ylabel('\rho','fontsize',16)
set(gca,'fontsize',16)
ylim([0,0.75])
saveas(gcf,['~/Desktop/',strrep(strrep(odes{ode_num},'.mat',''),'wendydata_',''),'_res.png'])

end