%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function V = get_VVp_svd(mt,t,K,phifun,center_scheme)
    dt = mean(diff(t));
    M = length(t);
%     v = phifun(linspace(-1,1,2*mt+1))*dt;
    Cfs = phi_weights(phifun,mt,1);
    v = Cfs(1,:)*dt;

    if isequal(center_scheme,'uni')
      gap = max(1,floor((M-2*mt)/K));
      %V = convmtx(v,M-2*mt);
      %V = V(1:gap:end,:);
      %Vp = convmtx(vp,M-2*mt);
      %Vp = V(1:gap:end,:);
      diags = 0:gap:M-2*mt-1;
      diags = diags(1:min(K,end));
        
      V = zeros(length(diags),M);
      for j=1:length(diags)
            V(j,gap*(j-1)+1:gap*(j-1)+2*mt+1) = v;
       end
    elseif isequal(center_scheme,'random')
        gaps = randperm(M-2*mt,K);        
        V = zeros(K,M);
        for j=1:K
            V(j,gaps(j):gaps(j)+2*mt) = v;
        end
    elseif isequal(class(center_scheme),'double')
        center_scheme = unique(max(min(center_scheme,M-mt),mt+1));
        K = length(center_scheme);
        V = zeros(K,M);
        for j=1:K
            V(j,center_scheme(j)-mt:center_scheme(j)+mt) = v;
        end        
    end
end


