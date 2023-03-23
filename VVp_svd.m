%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function [V,Vp] = VVp_svd(V,K_min,t,toggle_VVp_svd)
    m = length(t);
    dt = mean(diff(t));
    [U,S,~] = svd(V',0);
    sings = diag(S);
    if toggle_VVp_svd>0
        s = find(cumsum(sings.^2)/sum(sings.^2)>toggle_VVp_svd^2,1);
        if isempty(s)
            s = min(K,size(V,1));
        end
    else
        corner_data = cumsum(sings)/sum(sings);
        s = getcorner(corner_data,(1:length(corner_data))');%
        s = min(max(K_min,s),size(V,1));
%         plot(1:length(corner_data),corner_data,s,corner_data(s),'o'); drawnow
    end
    inds = 1:s;
    V = U(:,inds)'*dt;

    Vp = V';
    Vp_hat = fft(Vp);
    if mod(m,2)==0
        k = [0:m/2 -m/2+1:-1]';
    else
        k = [0:floor(m/2) -floor(m/2):-1]';
    end
    Vp_hat = ((2*pi/m/dt)*1i*k).*Vp_hat;
    % For odd derivatives there is a loss of symmetry
    if mod(m,2)==0
        Vp_hat(m/2) = 0;        
    end
    Vp = real(ifft(Vp_hat))';

%     M = triu(abs(V*Vp'+Vp*V')/dt)';
%     M(logical(speye(size(M,1)))) = 0;
%     max_ind = find(max(M,[],2) > 10^-14,1);
% %     max_ind = find(vecnorm(M,2,2) > 10^-14,1);
%     if isempty(max_ind)
%         max_ind = size(M,1);
%     else 
%         max_ind = max(4,max_ind);
%     end
%     inds_final = 1:max_ind;
% %     inds_final = abs(V*Vp'+Vp*V')/dt < 10^-10;
%     V = V(inds_final,:);
%     Vp = Vp(inds_final,:);
%     V = V./vecnorm(Vp,2,2)*sqrt(dt);
%     Vp = Vp./vecnorm(Vp,2,2)*sqrt(dt);

end

%     s = findchangepts(log10(diag(S)));% 
%     figure(15);clf; plot(log10(diag(S)/S(1,1))); hold on; plot(s,log10(S(s,s)/S(1,1)),'o')
%     s = find(diag(S)/S(1,1)<10^-1,1);
%     inds = randperm(size(V,1),min(K,size(V,1)));


function tstarind = getcorner(Ufft,xx)
    NN = length(Ufft);
    Ufft = Ufft/max(abs(Ufft))*NN;
    errs = zeros(NN,1);
    for k=1:NN
        [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k);
%         plot(1:k,L1,k:length(xx),L2,1:length(xx),Ufft,'b-o'); drawnow;
%         errs(k) = sqrt(sum(((L1-Ufft_av1)./Ufft_av1).^2) + sum(((L2-Ufft_av2)./Ufft_av2).^2)); % relative l2   
        errs(k) = sum(abs((L1-Ufft_av1)./Ufft_av1)) + sum(abs((L2-Ufft_av2)./Ufft_av2)); % relative l1   
%         errs(k) = m1^2/(1+m1^2)*norm(1:k)^2 + 1/(1+m1^2)*norm(Ufft_av1-b1)^2 + m2^2/(1+m2^2)*norm(k:NN)^2 + 1/(1+m2^2)*norm(Ufft_av2-b2)^2; % relative l2 
%         errs(k) = norm(L1-Ufft_av1,2)/norm(Ufft_av1,2)+norm(L2-Ufft_av2,2)/norm(Ufft_av2,2); % relative l2 
%         errs(k) = norm(L1-Ufft_av1,inf)/norm(Ufft_av1,inf)+norm(L2-Ufft_av2,inf)/norm(Ufft_av2,inf); % relative l2 
%         A = [ones(NN+1,1) [-m1*ones(length(Ufft_av1),1);-m2*ones(length(Ufft_av2),1)]]; % total least squares v2
%         c = [-b1*ones(length(Ufft_av1),1);-b2*ones(length(Ufft_av2),1)];
%         z = [[Ufft_av1(:) (1:k)']*(m1-1)/(m1^2+1); [Ufft_av2(:) (k:NN)']*(m2-1)/(m2^2+1)];
%         errs(k) = norm(sum(A.*z,2)-c);        
        
    end
    [~,tstarind] = min(errs);
end

function [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k)

   NN = length(Ufft);
   subinds1 = 1:k;
   subinds2 = k:NN;
   Ufft_av1 = Ufft(subinds1);
   Ufft_av2 = Ufft(subinds2);
   
   [m1,b1,L1]=lin_regress(Ufft_av1,xx(subinds1));
   [m2,b2,L2]=lin_regress(Ufft_av2,xx(subinds2));

end

function [m,b,L]=lin_regress(U,x)

   m = (U(end)-U(1))/(x(end)-x(1));
   b = U(1)-m*x(1);
   L = U(1)+m*(x-x(1));

%    A = [0*x(:)+1 x(:)];
%    c = A \ U(:);
%    b = c(1); m = c(2); L = A*c;

end