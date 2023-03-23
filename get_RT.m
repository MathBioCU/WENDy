%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function [RT,L0,Cov,diag_reg] = get_RT(L0,L1,w,diag_reg)
    dims = size(L1);
    if ~all(~w)
        L0 = L0 + reshape(reshape(permute(L1,[3 1 2]),dims(3),[]).'*w,dims(1),[]);
    end
    Cov = L0*(L0');
    RT = chol((1-diag_reg)*Cov+diag_reg*diag(diag(Cov)));
    RT = RT';
end