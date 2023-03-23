%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function sig = estimate_sigma(f)

    k = 6;
    C = fdcoeffF(k,0,-k-2:k+2);
    filter = C(:,end);
    filter = filter/norm(filter,2);
    sig = rms(reshape(conv2(filter(:),1,f,'valid'),[],1));

end