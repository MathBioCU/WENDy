function [V,Vp] = get_VVp(mt,t,max_d,K,phifun,center_scheme)
    dt = mean(diff(t));
    M = length(t);
    if ~isempty(phifun)
        Cfs = phi_weights(phifun,mt,max_d);
    else
        Cfs = [[0 0];fdcoeffF(1, t(mt), t(1:mt*2-1));[0 0]]';
        Cfs(2,:) = -Cfs(2,:)*(mt*dt);
    end
    v = Cfs(end-1,:)*(mt*dt)^(-max_d+1)*dt;
    vp = Cfs(end,:)*(mt*dt)^(-max_d)*dt;

    if isequal(center_scheme,'uni')
      gap = max(1,floor((M-2*mt)/K));
      %V = convmtx(v,M-2*mt);
      %V = V(1:gap:end,:);
      %Vp = convmtx(vp,M-2*mt);
      %Vp = V(1:gap:end,:);
      diags = 0:gap:M-2*mt-1;
      diags = diags(1:min(K,end));
        
      V = zeros(length(diags),M);
      Vp = zeros(length(diags),M);
      for j=1:length(diags)
            V(j,gap*(j-1)+1:gap*(j-1)+2*mt+1) = v;
            Vp(j,gap*(j-1)+1:gap*(j-1)+2*mt+1) = vp;
       end
    elseif isequal(center_scheme,'random')
        gaps = randperm(M-2*mt,K);        
        V = zeros(K,M);
        Vp = zeros(K,M);
        for j=1:K
            V(j,gaps(j):gaps(j)+2*mt) = v;
            Vp(j,gaps(j):gaps(j)+2*mt) = vp;
        end
    elseif isequal(class(center_scheme),'double')
        center_scheme = unique(max(min(center_scheme,M-mt),mt+1));
        K = length(center_scheme);
        V = zeros(K,M);
        Vp = zeros(K,M);
        for j=1:K
            V(j,center_scheme(j)-mt:center_scheme(j)+mt) = v;
            Vp(j,center_scheme(j)-mt:center_scheme(j)+mt) = vp;
        end        
    end

%     figure(1)
%     n=1;
%     for i=1:4
%         subplot(2,2,i)
%         plot(Vp_cell{1}(n+i,:),'.-')
%     end

end


