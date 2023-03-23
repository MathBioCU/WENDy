%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function Jac = get_jac(Jac_mat,w,features)
    param_length_vec = cellfun(@(x)length(x),features);
    w_cell = vec2cell(w,param_length_vec);
    Jac = pagemtimes(blkdiag(w_cell{:})',Jac_mat);
end