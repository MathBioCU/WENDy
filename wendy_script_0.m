%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

%% generate data

[t,x,features,params,x0,true_vec,rhs_p] = load_ode();
% [t,x,features,~,~,~,~] = load_ode(); % see load_ode for structurs of t,x,features

%% subsample timepoints

subsamp = 1;                       % subsample data in time
tobs = t(1:subsamp:end); 
xsub = x(1:subsamp:end,:);
[M,nstates] = size(xsub);

%% add noise

rng('shuffle');
noise_ratio = 0.1;
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
K_min = length([features{:}]);
mt_max = max(floor((M-1)/2)-K_min,1);
mt_min = rad_select(tobs,xobs,phifun,1,submt,0,1,2,mt_max,[]);
mt_cell = cellfun(@(x,y) [x,{y}], repmat({{phifun,meth}},length(mt_params),1),num2cell(mt_params(:)),'uni',0);

%% post-processing options

toggle_plot = 1;
toggle_ddd = 1;

%% run WENDy

use_true = 0;
tic;
[w_hat,res,res_0,w_hat_its,V_cell,Vp_cell,...
    Theta_cell,mt,xobs,Jac_mat,G_0,b_0,RT,stdW,mseW,CovW] = wendy_fcn_0(...
    xobs,tobs,features,toggle_smooth,...
    mt_cell,mt_min,mt_max,K_min,K_max,center_scheme,toggle_VVp_svd,...
    w0,optim_params,iter_diff_tol,max_iter,diag_reg,pvalmin,check_pval_it);
total_time = toc;

%% display results

if exist('true_vec','var')
    err_norm = 2;
    errs = arrayfunvec(w_hat_its,@(w)norm(w-true_vec,err_norm)/norm(true_vec,err_norm),1);
    err_wendy = errs(end);    
    disp(['err WENDy:',num2str(err_wendy)])
    disp('parameters: true, WENDy')
    disp([true_vec w_hat])
    G = RT \ G_0;
    b = RT \ b_0;
    res_true = G*true_vec-b;
    res_0_true = G_0*true_vec-b_0;
    err_NLS = 0;
    total_time_nls = 0;
    param_nls = [];
    display_wendy_results;
else
    display_wendy_results_0;
end

%%% define ODE
function [t,x,features,params,x0,true_vec,rhs_p] = load_ode()
    %%% SIRS
    numeq = 3;
    beta = 0.5; gamma = 0.1; delta = 0.1;
    x0 = [0.99 0.01 0]';
    t = 0:0.25:80;
    
    features = cell(numeq,1);
    features{1} = {@(S,I,R) S.*I, @(S,I,R) R};
    features{2} = {@(S,I,R) S.*I, @(S,I,R) I};
    features{3} = {@(S,I,R) I,@(S,I,R) R};
    params = {[-beta delta],[beta -gamma],[gamma -delta]};
    
    tol_ode = 10^-10;
    rhs_p = @(x,params) rhs_fun(features,params,x);
    true_vec = [params{:}]';
    options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
    [t,x]=ode45(@(t,x)rhs_p(x,params),t,x0,options_ode_sim);
end