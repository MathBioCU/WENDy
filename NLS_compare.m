function out=NLS_compare(xobs,params,rhs_pv,tobs,x0_nls,tol_ode_nls,odemethod_nls,thresh)
    yobs = get_sol(params,rhs_pv,tobs,x0_nls,tol_ode_nls,odemethod_nls,thresh);
    out = norm(vecnorm(xobs(1:size(yobs,1),:)-yobs))^2;
end

function x = get_sol(params,rhs,tspan,x0,tol_ode,odemethod,thresh)
    options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)),'Events',@(T,Y)myEvent(T,Y,thresh));
    if isequal(odemethod,'ode15s')
        [t,x]=ode15s(@(t,x)rhs(x,params),tspan,x0,options);  
    elseif or(isempty(odemethod),isequal(odemethod,'ode45'))
        [t,x]=ode45(@(t,x)rhs(x,params),tspan,x0,options);
    end
end

function [value, isterminal, direction] = myEvent(T, Y,thresh)
    value      = norm(Y(:),inf) >= thresh;
    isterminal = 1;   % Stop the integration
    direction  = 0;
end