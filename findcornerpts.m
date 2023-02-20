function [corner,sig_est] = findcornerpts(xobs,t)
    t = t(:);
    T = length(t);
    wn = ((0:T-1)-floor(T/2))'*(2*pi)/range(t);
    xx = wn(1:ceil(end/2));
    NN = length(xx);
    Ufft = mean(abs(fftshift(fft(xobs))/sqrt(2*NN)),2);
    Ufft = Ufft(1:ceil(T/2),:);
    [~,Umax] = max(Ufft);
    tstarind1 = getcorner(cumsum(abs(Ufft(1:Umax))),xx(1:Umax));
%     tstarind2 = getcorner(log(abs(Ufft(1:Umax))),xx(1:Umax));
%     tstarind = floor((tstarind1+tstarind2)/2);
    tstarind = tstarind1;
%     tstarind = tstarind2;
    tstar = -xx(tstarind);
    corner = [tstar max(Umax-tstarind+1,1)];
    sig_est = rms(Ufft(1:min(corner(:,2)),:));
end