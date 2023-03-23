%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function [m,sigma_est,its] = get_optimal_SMAF(x,fx_obs,max_points,init_m_fac,max_filter_fac,expand_fac,maxits,deriv_tol,verbose)

    if isempty(max_points)
        max_points = 10^5;
    end
    if isempty(init_m_fac)
        init_m_fac = 200;
    end
    if isempty(max_filter_fac)
        max_filter_fac = 8;
    end
    if isempty(expand_fac)
        expand_fac = 2;
    end
    if isempty(maxits)
        maxits = 100;
    end
    if isempty(deriv_tol)
        deriv_tol = 10^-6;
    end
    if isempty(verbose)
        verbose = 0;
    end

    sigma_est = estimate_sigma(fx_obs,6,1);
    subsamp = max(floor(length(x)/max_points),1);  
    
    fx_subsamp = fx_obs(1:subsamp:end,:);
    M = size(fx_subsamp,1);
    dx_subsamp = mean(diff(x(1:subsamp:end)));
    
    m = ceil(M/init_m_fac);
    max_filter_width=floor(M/max_filter_fac);
    
    its = 1; check=1; 
    m = min(m,max_filter_width);

    while and(check>0,its<maxits)
        if verbose
            tic
        end
        [~,A] = build_poly_kernel(2,@(x) x*0+1,min(max(floor(m*expand_fac),3),floor((M-1)/2)),dx_subsamp,0);
        if size(fx_subsamp,2)==1
            d = 2*mean(abs(conv(fx_subsamp,A(3,:),'valid')));
        else
            d = 2*mean(reshape(abs(conv2(A(3,:),1,fx_subsamp,'valid')),[],1));
        end
        C = sigma_est^2/((d+deriv_tol)^2*dx_subsamp^4/144);
        mnew = min( floor((fzero(@(n) n.^5-n.^3-C,1)-1)/2), max_filter_width);
    
        check = abs(m-mnew);
        m = mnew;
        its = its+1;
        if verbose
            disp([toc m d])
        end
    end
    m = m*subsamp;
    
end

function [f,A]=build_poly_kernel(deg,k,n,dx,max_dx)
    x = (-n:n)'*dx;
    X = x.^(0:deg);
    K = k(x/(n*dx));
    K = K/norm(K,1);
    A = pinv(sqrt(K).*X).*sqrt(K)';
    M = [diag(factorial(0:max_dx)) zeros(max_dx+1,deg-max_dx)];
    f = M*A;
end

function sig = estimate_sigma(f,k,dim)

    C = fdcoeffF(k,0,-k-2:k+2);
    filter = C(:,end);
    filter = filter/norm(filter,2);
    if dim>1
        sig = rms(reshape(convn(permute(f,[dim 1:dim-1 dim+1:length(size(f))]),filter(:),'valid'),[],1));
    else
        if size(f,2)==1
            sig = rms(conv(f,filter(:),'valid'));
        else
            sig = rms(reshape(conv2(filter(:),1,f,'valid'),[],1));
        end
    end

end

function C = fdcoeffF(k,xbar,x)

% Compute coefficients for finite difference approximation for the
% derivative of order k at xbar based on grid values at points in x.
%
% This function returns a row vector c of dimension 1 by n, where n=length(x),
% containing coefficients to approximate u^{(k)}(xbar), 
% the k'th derivative of u evaluated at xbar,  based on n values
% of u at x(1), x(2), ... x(n).  
%
% If U is a column vector containing u(x) at these n points, then 
% c*U will give the approximation to u^{(k)}(xbar).
%
% Note for k=0 this can be used to evaluate the interpolating polynomial 
% itself.
%
% Requires length(x) > k.  
% Usually the elements x(i) are monotonically increasing
% and x(1) <= xbar <= x(n), but neither condition is required.
% The x values need not be equally spaced but must be distinct.  
%
% This program should give the same results as fdcoeffV.m, but for large
% values of n is much more stable numerically.
%
% Based on the program "weights" in 
%   B. Fornberg, "Calculation of weights in finite difference formulas",
%   SIAM Review 40 (1998), pp. 685-691.
%
% Note: Forberg's algorithm can be used to simultaneously compute the
% coefficients for derivatives of order 0, 1, ..., m where m <= n-1.
% This gives a coefficient matrix C(1:n,1:m) whose k'th column gives
% the coefficients for the k'th derivative.
%
% In this version we set m=k and only compute the coefficients for
% derivatives of order up to order k, and then return only the k'th column
% of the resulting C matrix (converted to a row vector).  
% This routine is then compatible with fdcoeffV.   
% It can be easily modified to return the whole array if desired.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)


n = length(x);
if k >= n
   error('*** length(x) must be larger than k')
   end

m = k;   % change to m=n-1 if you want to compute coefficients for all
         % possible derivatives.  Then modify to output all of C.
c1 = 1;
c4 = x(1) - xbar;
C = zeros(n-1,m+1);
C(1,1) = 1;
for i=1:n-1
  i1 = i+1;
  mn = min(i,m);
  c2 = 1;
  c5 = c4;
  c4 = x(i1) - xbar;
  for j=0:i-1
    j1 = j+1;
    c3 = x(i1) - x(j1);
    c2 = c2*c3;
    if j==i-1
      for s=mn:-1:1
        s1 = s+1;
        C(i1,s1) = c1*(s*C(i1-1,s1-1) - c5*C(i1-1,s1))/c2;
        end
      C(i1,1) = -c1*c5*C(i1-1,1)/c2;
      end
    for s=mn:-1:1
      s1 = s+1;
      C(j1,s1) = (c4*C(j1,s1) - s*C(j1,s1-1))/c3;
      end
    C(j1,1) = c4*C(j1,1)/c3;
    end
  c1 = c2;
end            % last column of c gives desired row vector
end
