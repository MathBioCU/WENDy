%%% lorenz
% features = cell(3,1);
% features{1} = {@(x,y,z) y, @(x,y,z) x};
% features{2} = {@(x,y,z) x, @(x,y,z) x.*z, @(x,y,z) y};
% features{3} = {@(x,y,z) x.*y, @(x,y,z) z};
% params = {[10 -10],[28 -1 -1],[1 -8/3]};
% x0 = [-8 7 27];
% t = [0:0.02:12];

%% SIR

b = 0.003;
g = 0.25;
features = cell(3,1);
features{1} = {@(S,I,R) S.*I};
features{2} = {@(S,I,R) S.*I, @(S,I,R) I};
features{3} = {@(S,I,R) I};
params = {-b,[b -g],g};
x0 = [499 1 0];
t = 0:0.04:40;

%% SIRS

beta = 0.5; gamma = 0.1; delta = 0.1;
x0 = [0.99 0.01 0]';
t = 0:0.1:100;
features = cell(3,1);
features{1} = {@(S,I,R) S.*I, @(S,I,R) R};
features{2} = {@(S,I,R) S.*I, @(S,I,R) I};
features{3} = {@(S,I,R) I,@(S,I,R) R};
params = {[-beta delta],[beta -gamma],[gamma -delta]};

%% features, params, x0, t ---->  data

tic;
tol_ode = 10^-8;
rhs_p = @(x,params) rhs_fun(features,params,x);
true_vec = [params{:}]';
options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
[t,x]=ode45(@(t,x)rhs_p(x,params),t,x0,options_ode_sim);
disp(['sim time=',num2str(toc)])
plot(t,x)

% dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/'; 
% filename = 'gyroceptron_r_ep05.mat';
% load([dr,filename]);
