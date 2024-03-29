%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function Jac_mat = build_Jac_sym(features,xobs)
    [M,nstates] = size(xobs);
    features = [features{:}];
    J = length(features);
    Jac_mat = zeros(J,nstates,M);
    xobs_cell = mat2cell(num2cell(xobs'),nstates,ones(1,M));
    X = str2sym(strcat('x',num2str((1:nstates)')));
    Xc = sym2cell(X);
    for j=1:J
        f = features{j};
        g = gradient(f(Xc{:}),X);
        G = matlabFunction(g,'vars',Xc);
        z = cell2mat(cellfun(@(x)G(x{:}),xobs_cell,'uni',0));
        Jac_mat(j,:,:) = z;
    end
end