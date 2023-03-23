%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

%%% general radius finder for tf phifun on [-1,1]. Each 
%%% meth must depend on only one positive parameter p
%%% xobs must be a column vector, so acts on each component of xobs

function mt = get_rad(xobs,tobs,phifun,meth,p,mt_min,mt_max)
    if any([isempty(phifun) isempty(meth) isempty(p)])
        mt = 0;
    else
        if isequal(meth,'direct')
            mt = p;
        elseif isequal(meth,'FFT')
            [mt,~,~,~] = findcorners(xobs,tobs,[],p,phifun);
        elseif isequal(meth,'timefrac')
            mt = floor(length(tobs)*p);
        elseif isequal(meth,'mtmin')
            mt = p*mt_min;
        end
        mt = min(mt,mt_max);
    end
end