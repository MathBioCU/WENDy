% t = [0:0.025:30]; %vanderpol (up to 20% noise)
% t = [0:0.08:24]; %duffing (up to 20% noise)
% t = [0:0.02:6]; %lotka-volterra (up to 20% noise)
% t = [0:0.02:6]; %lorenz (up to 20% noise)
% t = [0:0.08:24]; %rossler (up to 15% noise)
% t = [0:0.05:20]; %oregonator (up to 2% noise)
% t = [0:0.005:10]; %Hindmarsh-Rose (up to 4% noise)
% t = [0:0.15:30]; %pendulum (up to 20% noise)
% t = [0:0.05:15]; %cubicOsc (up to 20% noise)
% t = [0:0.04:4]; %Gompertz (up to 10% noise, but jac correct makes it worse)
% t = [0:0.1:20]; %FitzHugh (up to 15% noise)
% t = [0:0.1:500]; %gyroceptron 

function [features,params,x0,t] = get_ode_features(ode_name)
    if isequal(ode_name,'Logistic_Growth')
        features = {{@(x) x, @(x) x.^2}};
        params = {[1 -1]};
        x0 = [0.01];
        t = [0:0.04:10];
    elseif isequal(ode_name,'Gompertz')
        features = {{@(x) x, @(x) x.*log(abs(x))}};
        params = {[2 -2]};
        x0 = [0.01];
        t = [0:0.02:5];
    elseif isequal(ode_name,'rational')
        features = {{@(x) 1./(1+x.^2), @(x) x./(1+x.^2)}};
        params = {[1 -1]};
        x0 = [-5];
        t = [0:0.04:4];
    elseif isequal(ode_name,'Linear')
        features = cell(2,1);
        features{1} = {@(x,y) x,@(x,y) y};
        features{2} = {@(x,y) x,@(x,y) y};
        params = {[-0.1 2],[-2 -0.1]};
        x0 = [0 2];
        t = [0:0.05:20];
    elseif isequal(ode_name,'Van_der_Pol')
        features = cell(2,1);
        features{1} = {@(x,y) y};
        features{2} = {@(x,y) x,@(x,y) y,@(x,y) x.^2.*y};
        params = {1,[-1 4 -4]};
        x0 = [1 1];
        t = [0:0.08:20];
    elseif isequal(ode_name,'Duffing')
        features = cell(2,1);
        features{1} = {@(x,y) y};
        features{2} = {@(x,y) x, @(x,y) y,@(x,y) x.^3};
        params = {1,[-0.2 -0.2 -1]};
        x0 = [1 1];
        t = [0:0.08:24];
    elseif isequal(ode_name,'Lotka_Volterra')
        features = cell(2,1);
        features{1} = {@(x,y) x, @(x,y) x.*y};
        features{2} = {@(x,y) y, @(x,y) x.*y};
        params = {[3 -1],[-6 1]};
        x0 = [1 1];    
        t = [0:0.02:5];
    elseif isequal(ode_name,'cubicOsc')
        features = cell(2,1);
        features{1} = {@(x,y) x.^3, @(x,y) y.^3};
        features{2} = {@(x,y) x.^3, @(x,y) y.^3};
        params = {[-0.1 2],[-2 -0.1]};
        x0 = [0 2]; 
        t = [0:0.05:15];
    elseif isequal(ode_name,'FitzHugh-Nagumo')
        features = cell(2,1);
        features{1} = {@(x,y) x, @(x,y) x.^3, @(x,y) y};
        features{2} = {@(x,y) x, @(x,y) x*0+1, @(x,y) y};
        params = {[3 -3 3],[-1/3 0.34/3 0.2/3]};
        x0 = [0 0.1]; 
        t = linspace(0,25,401);
    elseif isequal(ode_name,'pendulum')
        features = cell(2,1);
        features{1} = {@(x,y) y};
        features{2} = {@(x,y) sin(x)};
        params = {1,-1};
        x0 = [pi-1*pi/16 0];
        t = [0:0.15:30];
    elseif isequal(ode_name,'Lorenz')
        features = cell(3,1);
        features{1} = {@(x,y,z) y, @(x,y,z) x};
        features{2} = {@(x,y,z) x, @(x,y,z) x.*z, @(x,y,z) y};
        features{3} = {@(x,y,z) x.*y, @(x,y,z) z};
        params = {[10 -10],[28 -1 -1],[1 -8/3]};
        x0 = [-8 7 27];
        t = [0:0.02:6];
    elseif isequal(ode_name,'Rossler')
        features = cell(3,1);
        features{1} = {@(x,y,z) y, @(x,y,z) z};
        features{2} = {@(x,y,z) x, @(x,y,z) y};
        features{3} = {@(x,y,z) x*0+1, @(x,y,z) x.*z, @(x,y,z) z};
        params = {[-1 -1],[1 0.2],[0.2 1 -5.7]};
        x0 = [3 5 0];
        t = [0:0.08:20];
    elseif isequal(ode_name,'Oregonator')
        features = cell(3,1);
        features{1} = {@(x,y,z) y, @(x,y,z) x.*y, @(x,y,z) x, @(x,y,z) x.^2};
        features{2} = {@(x,y,z) y, @(x,y,z) x.*y, @(x,y,z) z};
        features{3} = {@(x,y,z) x, @(x,y,z) z};
        params = {[0.5 -5 5 -1],[-0.5 -5 0.2375],[10 -0.5]};
        x0 = [1 0.2 50];
        t = [0:0.05:13];
    elseif isequal(ode_name,'Hindmarsh-Rose')
        features = cell(3,1);
        features{1} = {@(x,y,z) y, @(x,y,z) x.^3, @(x,y,z) x.^2, @(x,y,z) z};
        features{2} = {@(x,y,z) x*0+1, @(x,y,z) x.^2, @(x,y,z) y};
        features{3} = {@(x,y,z) x, @(x,y,z) x*0+1, @(x,y,z) z};
        params = {[10 -10 30 -10],[10 -50 -10],[0.04 0.0319 -0.01]};
        x0 = [-1.3095   -7.5901   -0.2020];
        t = [0:0.0025:10];
    elseif isequal(ode_name,'gyroceptron')
        ep= 0.05;
        features = cell(4,1);
        features{1} = {@(q1,p1,q2,p2) p1};
        features{2} = {@(q1,p1,q2,p2) q1,@(q1,p1,q2,p2) q2.*sin(2*q1+2*q2)+2*q1.*q2.*cos(2*q1+2*q2)};
        features{3} = {@(q1,p1,q2,p2) p2};
        features{4} = {@(q1,p1,q2,p2) q2,@(q1,p1,q2,p2) q1.*sin(2*q1+2*q2)+2*q1.*q2.*cos(2*q1+2*q2)};
        params = {[1],[-1 ep],[ep],[-ep -ep]};
        x0 = [0 1 0 1];
        t = [0:0.1:500];
    elseif isequal(ode_name,'gyroceptron_r')
        ep = 0.05;
        features = cell(2,1);
        features{1} = {@(q2,p2) p2};
        features{2} = {@(q2,p2) q2, @(q2,p2) cos(2*q2), @(q2,p2) q2.*sin(2*q2)};
        mu = 0.5;
        params = {[ep],[-ep -ep*sqrt(2*mu)*besselj(1,2*sqrt(2*mu)) ep*2*sqrt(2*mu)*besselj(1,2*sqrt(2*mu))]};
        x0 = [0 1 0 1];
        t = [0:0.1:500];
    elseif isequal(ode_name,'biochemM1')
        % Bayesian ranking of biochemical system models
        features = cell(5,1);
        Km = 0.3;
        features{1} = {@(S,dS,R,RS,Rpp) S, @(S,dS,R,RS,Rpp) S.*R, @(S,dS,R,RS,Rpp) RS};
        features{2} = {@(S,dS,R,RS,Rpp) S};
        features{3} = {@(S,dS,R,RS,Rpp) S.*R, @(S,dS,R,RS,Rpp) RS, @(S,dS,R,RS,Rpp) Rpp./(Km + Rpp)};
        features{4} = {@(S,dS,R,RS,Rpp) S.*R, @(S,dS,R,RS,Rpp) RS};
        features{5} = {@(S,dS,R,RS,Rpp) RS, @(S,dS,R,RS,Rpp) Rpp./(Km + Rpp)};
        params = {[-0.07 -0.6 0.35],[0.07],[-0.6 0.05 0.017],[0.6 -0.35],[0.3 -0.017]};
        x0 = [1 0 1 0 1];
        t = [0:0.1:25];
    elseif isequal(ode_name,'alphapinene')
        features = cell(5,1);
        features{1} = {@(x1,x2,x3,x4,x5) x1};
        features{2} = {@(x1,x2,x3,x4,x5) x1};
        features{3} = {@(x1,x2,x3,x4,x5) x1, @(x1,x2,x3,x4,x5) x3, @(x1,x2,x3,x4,x5) x5};
        features{4} = {@(x1,x2,x3,x4,x5) x3};
        features{5} = {@(x1,x2,x3,x4,x5) x3, @(x1,x2,x3,x4,x5) x5};
        params = {400*(-5.93e-5-2.96e-5), 400*(5.93e-5), 400*([2.96e-5 -(2.05e-5+27.5e-5) 4e-5]), 400*(2.05e-5),400*([27.5e-5 -4e-5])};
        x0 = [100 0 0 0 0];
        t = [0:0.25:10^2];
    elseif isequal(ode_name,'lorenz96')
        features = cell(6,1);
        features{1} = {@(x1,x2,x3,x4,x5,x6) x2.*x6, @(x1,x2,x3,x4,x5,x6)x5.*x6, @(x1,x2,x3,x4,x5,x6) x1, @(x1,x2,x3,x4,x5,x6) x1*0+1};
        features{2} = {@(x1,x2,x3,x4,x5,x6) x3.*x1, @(x1,x2,x3,x4,x5,x6)x6.*x1, @(x1,x2,x3,x4,x5,x6) x2, @(x1,x2,x3,x4,x5,x6) x1*0+1};
        features{3} = {@(x1,x2,x3,x4,x5,x6) x4.*x2, @(x1,x2,x3,x4,x5,x6)x1.*x2, @(x1,x2,x3,x4,x5,x6) x3, @(x1,x2,x3,x4,x5,x6) x1*0+1};
        features{4} = {@(x1,x2,x3,x4,x5,x6) x5.*x3, @(x1,x2,x3,x4,x5,x6)x2.*x3, @(x1,x2,x3,x4,x5,x6) x4, @(x1,x2,x3,x4,x5,x6) x1*0+1};
        features{5} = {@(x1,x2,x3,x4,x5,x6) x6.*x4, @(x1,x2,x3,x4,x5,x6)x3.*x4, @(x1,x2,x3,x4,x5,x6) x5, @(x1,x2,x3,x4,x5,x6) x1*0+1};
        features{6} = {@(x1,x2,x3,x4,x5,x6) x1.*x5, @(x1,x2,x3,x4,x5,x6)x4.*x5, @(x1,x2,x3,x4,x5,x6) x6, @(x1,x2,x3,x4,x5,x6) x1*0+1};
        params = repmat({[1 -1 -1 8]},6,1);
        x0 = 8*ones(6,1); x0(1)=x0(1)+0.01;
        t = [0:0.02:10];
    else
        features={};
    end
end