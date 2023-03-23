%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function x0=x0est(x,x0_true,sigma_x0_nls)
    if isempty(x0_true)      
        diffx = [x(2,:)-x(1,:);(x(3:end,:)-x(1:end-2,:))/2];
        x0 = mean(x) -mean(cumsum(movmean(diffx,2,'endpoints','discard')));
    else
        x0 = x0_true+ randn(size(x0_true))*sigma_x0_nls.*abs(x0_true);
    end
end