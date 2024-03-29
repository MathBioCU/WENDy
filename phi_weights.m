%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger
 
function Cfs = phi_weights(phifun,m,maxd)
    xf = linspace(-1,1,2*m+1);
    x = xf(2:end-1);
    Cfs = zeros(maxd+1,2*m+1);
    syms y;
    f = @(y)phifun(y);
    for j=1:maxd+1
        Df = matlabFunction(diff(f(y),j-1));
        Cfs(j,2:end-1) = fillmissing(Df(x),'constant',Df(eps));
        inds = find(isinf(abs(Cfs(j,:))));
        for k=1:length(inds)
            Cfs(j,inds(k)) = Df(xf(inds(k))-sign(xf(inds(k)))*eps);
        end
    end
    Cfs = Cfs/norm(Cfs(1,:),2);
end
