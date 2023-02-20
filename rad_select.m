% inc = 1;
% sub = 3;
% q = 0;
% s=2;

function mt = rad_select(t0,y,phifun,inc,sub,q,s,m_min,m_max,pow)
    
    if isempty(phifun)
        mt = m_min;
    else
        [M,nstates] = size(y);
        dt = mean(diff(t0));
        t = 0:dt/inc:(M-1)*dt;
        
        if q>0
        %     ks = -inc*floor(M/2):(ceil(M/2)*inc-1);
        %     prox_u = @(k) 1/sqrt(M*dt)./(1+0.1*abs(k).^q);
        %     prox_u_vec = prox_u(ks);
            prox_u = @(t) exp(-abs(t-t(floor(end/2))).^q);
            prox_u_vec =  dt/inc/sqrt(M*dt)*fftshift(fft(prox_u(t)));
        else
            if inc>1
                prox_u_vec =  dt/inc/sqrt(M*dt)*fftshift(fft(interp1(t0,y,t,'spline')));
            elseif inc ==1
                prox_u_vec =  dt/sqrt(M*dt)*fftshift(fft(y));
            end
        end
        errs = [];
        ms = [];
        
        for m=m_min:m_max
            t_phi = linspace(-1+dt/inc,1-dt/inc,2*inc*m-1);
            Qs = 1:floor(s*inc*m):length(t)-2*inc*m;
            errs_temp = zeros(nstates,length(Qs));
            for Q=1:length(Qs)
                phi_vec = zeros(1,length(t));
                phi_vec(Qs(Q):Qs(Q)+length(t_phi)-1)=phifun(t_phi);
                phi_vec = phi_vec/norm(phi_vec,2);
                for nn=1:nstates
                    phiu_fft = (dt/sqrt(M*dt))*fft(phi_vec(:).*y(:,nn));
                    alias = phiu_fft(1:floor(M/sub):floor(inc*M/2));
                    errs_temp(nn,Q) = 2*(2*pi/sqrt(M*dt))*sum((0:length(alias)-1).*imag(alias(:)'));
                end
            end
            check1 = rms(errs_temp(:));
            errs = [errs check1];
            ms = [ms m];
        end
        
        if isequal(class(pow),'char')
            b = findchangepts(log(errs),'statistic',pow);
            if isempty(b)
                [~,b]=min(errs.*ms.^(1/2));
            else
                b=b(1);
            end
        elseif isempty(pow)
            b = getcorner(log(errs),ms);
        else
            [~,b]=min(errs.*ms.^pow);
        end
        mt = ms(b);
    %     figure(15);
    %     semilogy(ms,errs,'r.-',ms(b),errs(b),'bx')
    %     title(ms(b))
    %     ylim([min(errs)/10 max(errs)*10])
    end
end


%             phiu_fft = 1/sqrt(M*dt)*conv(phifft,prox_u_vec(:,nn),'same');    
%                 if ~mod(M/2,2)
%                     alias1 = phiu_fft(M/2+1:-floor(M/sub):1);
%                     alias_k1 = 0:-1:-length(alias1)+1;
%                     alias2 = phiu_fft(M/2+1+floor(M/sub):floor(M/sub):end);
%                     alias_k2 = 1:length(alias2);
%                     alias = [alias1(:);alias2(:)];
%                     alias_k = [alias_k1(:);alias_k2(:)];
%                 elseif mod(M/2,2)
%                     alias1 = phiu_fft((M-1)/2+1:-floor(M/sub):1);
%                     alias_k1 = 0:-1:-length(alias1)+1;
%                     alias2 = phiu_fft((M-1)/2+1+floor(M/sub):floor(M/sub):end);
%                     alias_k2 = 1:length(alias2);
%                     alias = [alias1(:);alias2(:)];
%                     alias_k = [alias_k1(:);alias_k2(:)];
%                 end     
%                 errs_temp(nn,Q) = 2*pi*1i/sqrt(M*dt)*sum(alias.*alias_k);
% 

%     t_phi = linspace(-1+dt/inc,1-dt/inc,2*inc*m-1);
%     Qs = inc*m+1:floor(s*inc*m):length(t)-inc*m;
%     errs_temp = zeros(nstates,length(Qs));
%     for Q=1:length(Qs)
%         phi_vec = [zeros(1,floor(length(t)/2)-inc*m-Qs(Q)) phifun(t_phi) zeros(1,ceil(length(t)/2)-inc*m+1+Qs(Q))];
%         phi_vec = phi_vec/norm(phi_vec,1);
%         phifft = dt/inc/sqrt(M*dt)*fftshift(fft(phi_vec));            
%         for nn=1:nstates
%             phiu_fft = 1/sqrt(M*dt)*conv(phifft,prox_u_vec(:,nn),'same');    
%             alias = ifftshift(phiu_fft);
%             alias = alias(1:floor(M/sub):floor(inc*M/2));
%             errs_temp(nn,Q) = 2*2*pi/sqrt(M*dt)*sum((0:length(alias)-1).*imag(alias(:)'));
%         end
%     end
