%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function [mts,pts,sig_ests,corners] = findcorners(xobs,t,tau,tauhat,phi_class)
    T = length(t);
    n = size(xobs,2);
    corners = zeros(n,2);
    sig_ests = zeros(n,1);
    mts = zeros(n,1);
    pts = zeros(n,1);
    
    if length(tauhat)==1
        tauhat = repmat(tauhat,1,n);
    end
    
    for nn= 1:n
        [corner,sig_est] = findcornerpts(xobs(:,nn),t);
        k = corner(2);
        if isequal(phi_class,1)
            l = @(m,k,N) log((2*m-1)./m.^2).*(4*pi^2*k^2*m.^2-3*N^2*tauhat(nn)^2)-2*N^2*tauhat(nn)^2*log(tau);
            mnew = fzero(@(m)l(m,k,T), [1 2/sqrt(tau)]);
            if mnew>T/2-1
                mnew = T/2/k;
            end
            mts(nn) = min(floor((T-1)/2),ceil(mnew)); 
            pts(nn) = max(2,floor(log(tau)/log(1-(1-1/mts(nn))^2)));
        elseif isequal(phi_class,2)
            mnew = 1+T*tauhat(nn)/2/pi/k*sqrt(-2*log(tau));
            mts(nn) = min(floor((T-1)/2),ceil(mnew)); 
            pts(nn) = 2*pi*k/tauhat(nn)/T;
        elseif isequal(class(phi_class), 'function_handle')
            mts(nn) = get_tf_support(phi_class,T,tauhat(nn),k);
            pts(nn) = NaN;
        end
        corners(nn,:)=corner;
        sig_ests(nn) = sig_est;
    end
end