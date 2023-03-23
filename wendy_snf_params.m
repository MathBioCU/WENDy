%% set and forget (iteration params)

%%% smooth data
toggle_smooth = 0;

%%% set optimization params
optim_params = {'meth','LS'};

%%% set weak integration
eta = 9;
phifuns = {@(x) exp(-eta*(1-x.^2).^(-1)), @(x) (1-x.^2).^eta};
center_scheme = 'uni';
toggle_VVp_svd = NaN; % 0, no SVD reduction; in (0,1), truncates Frobenious norm; NaN, truncates SVD according to cornerpoint of cumulative sum of singular values
submt = 2.1;

%%% set jacobian correction params
max_iter = 100;
iter_diff_tol = 10^-6;
err_norm = 2;
diag_reg = 10^-10; %% arbitrary low value to avoid warnings
check_pval_it = 10;
pvalmin = 10^-4;
w0 = []; % [], weak-form OLS initial guess