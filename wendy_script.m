%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

%% generate data

odes={'wendydata_Logistic_Growth.mat',...
    'wendydata_Lotka_Volterra.mat',...
    'wendydata_FitzHugh-Nagumo.mat',...
    'wendydata_Hindmarsh-Rose.mat',...
    'wendydata_biochemM1.mat'};

ode_num = 5;                       % select ODE from list above
load(odes{ode_num},'t','x','features','params','x0','true_vec','rhs_p');

subsamp = 4;                       % subsample data in time
tobs = t(1:subsamp:end); xsub = x(1:subsamp:end,:);
[M,nstates] = size(xsub);

%% add noise

rng(1);
noise_ratio = 0.05;
noise_dist = 0;
noise_alg = 0;
rng_seed = rng().Seed; rng(rng_seed);
[xobs,noise,~,sigma] = gen_noise(xsub,noise_ratio,noise_dist,noise_alg);

%% set WENDy params

%%% set-and-forget params
wendy_snf_params;

%%% set weak integration
phifun = phifuns{1};                % defined in wendy_snf_params.m
meth = 'mtmin';                     % 'mtmin','FFT','direct','timefrac'
mt_params = 2.^(0:3);               % see get_rad.m
K_max = 5000;
K_min = length(true_vec)*2;
mt_max = max(floor((M-1)/2)-K_min,1);
mt_min = rad_select(tobs,xobs,phifun,1,submt,0,1,2,mt_max,[]);
mt_cell = cellfun(@(x,y) [x,{y}], repmat({{phifun,meth}},length(mt_params),1),num2cell(mt_params(:)),'uni',0);

%% post-processing options

toggle_plot = 1;
toggle_ddd = 1;
toggle_nls = 0;

%% run WENDy

tic;
[w_hat,res,res_true,res_0,res_0_true,w_hat_its,errs,V_cell,Vp_cell,...
    Theta_cell,mt,x_jac,Jac_mat,G_0,b_0,RT,stdW,mseW,CovW] = wendy_fcn(...
    xobs,tobs,features,true_vec,toggle_smooth,...
    mt_cell,mt_min,mt_max,K_min,K_max,center_scheme,toggle_VVp_svd,...
    w0,optim_params,err_norm,iter_diff_tol,max_iter,diag_reg,pvalmin,check_pval_it);
total_time = toc;
err_windy = errs(end);

%% compare with NLS

tol_ode_nls = 10^-6;
bnds = {[],[]};
odemethod_nls = 'ode15s'; % 'ode15s', 'ode45'
sigma_x0_nls = 0;
param_init_true = 1;
sigma_param_init_nls = 0.25;

x0_nls = x0est(xobs,x0,sigma_x0_nls);
thresh = 10*max(vecnorm(xobs,inf,2));
rhs_pv = @(x,params) rhs_p(x,vec2cell(params,cellfun(@(x)length(x),features)));
fun = @(xobs,params) NLS_compare(xobs,params,rhs_pv,tobs,x0_nls,tol_ode_nls,odemethod_nls,thresh);
if param_init_true==1
    % uniform / sign preserving random initial guess
    param_init = true_vec+2*sqrt(3)*(rand(length(true_vec),1)-0.5).*abs(true_vec)*sigma_param_init_nls;
elseif param_init_true == 2
    % normally distributed random initial guess
    param_init = randn(size(true_vec))*sigma_param_init_nls;
elseif param_init_true == 3
    % WENDy initial guess
    param_init = w_hat;
elseif param_init_true==4
    % weak-form OLS initial guess
    param_init = G_0 \ b_0;
end
options_nls = optimoptions(@lsqnonlin,...
    'Algorithm','Levenberg-Marquardt',...
    'Display','iter',...
    'MaxIterations',500,'MaxFunctionEvaluations',2000,...
    'OptimalityTolerance',1e-8,'StepTolerance',1e-8);
if toggle_nls==1
    tic,
    param_nls = lsqnonlin(@(params)fun(xobs,params),param_init,bnds{:},options_nls);
    err_NLS = norm(true_vec-param_nls,err_norm)/norm(true_vec,err_norm);
    total_time_nls = toc;
else
    param_nls = true_vec*0;
    err_NLS = [];
    total_time_nls = [];
end

%% display results

display_wendy_results;