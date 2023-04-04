features = cell(3,1);
features{1} = {@(x,y,z) y, @(x,y,z) x};
features{2} = {@(x,y,z) x, @(x,y,z) x.*z, @(x,y,z) y};
features{3} = {@(x,y,z) x.*y, @(x,y,z) z};
params = {[10 -10],[28 -1 -1],[1 -8/3]};
x0 = [-8 7 27];
t = [0:0.02:12];

tol_ode = 10^-8;

rhs_p = @(x,params) rhs_fun(features,params,x);
true_vec = [params{:}]';
options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));

tic;
[t,x]=ode45(@(t,x)rhs_p(x,params),t,x0,options_ode_sim);
disp(['sim time=',num2str(toc)])
plot(t,x)