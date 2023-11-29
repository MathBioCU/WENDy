disp(['-----------------'])
K = size(V_cell{1},1);
disp(['time WENDy: ',num2str(total_time), '(sec)'])
disp(['Num its: ', num2str(size(w_hat_its,2))])
disp(['K,d,M:',num2str(K),',',num2str(nstates),',',num2str(length(tobs))])
disp(['test function radii: (',num2str(unique(mt)'),')'])
disp([' '])

%%% get rhs_wendy
ind = inf; % choose which weights to generate from
rhs_p = @(x,params) rhs_fun(features,params,x);
w_jac_cell = vec2cell(w_hat_its(:,min(ind,end)),cellfun(@(x)length(x),features));
rhs_wendy = @(x) rhs_p(x,w_jac_cell);

%%% get data-driven dynamics
if toggle_ddd
    x0 = xsub(1,:);
    tol_ode = 1e-12;
    options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
    [t_learned,x_learned]=ode45(@(t,x)rhs_wendy(x),tobs,x0,options_ode_sim);
end

%%% compute FFTs
tau = 10^-10;
tauhat = 1;
[~,~,~,corners] = findcorners(xobs,tobs,tau,tauhat,phifun);
xfft = abs(fft(xobs)); xfft = xfft./max(xfft);
phifft = cell2mat(cellfun(@(V)abs(fft(full(V(1,:)))),Vp_cell,'uni',0)); 
phifft = phifft./max(phifft,[],2);

%%% compute conf intervals
xflip = [1:length(w_hat) length(w_hat):-1:1];
c = 0.05; % <(100)c chance of not containing true val
stdW = max(sqrt(diag(CovW)),eps);
conf_int = arrayfun(@(x)norminv(1 - c/2,0,x),stdW);
confbounds = [w_hat-conf_int;flipud(w_hat+conf_int)];

fignum=1;
figure(fignum)
fignum=fignum+1;
h1=plot(1:length(w_hat),w_hat,'ro');hold on
for j=1:length(w_hat)
    rectangle('position',[j-0.25 w_hat(j)-conf_int(j) 0.5 2*conf_int(j)]);
    line([j-0.25 j+0.25],[w_hat(j) w_hat(j)])
end
hold off
legend(h1,{'{\bf w}_{WENDy}'},'box','off','location','best')
title([num2str((1-c)*100),'% confidence bounds'],'fontsize',9)
set(gca,'Xtick',1:length(w_hat))
xlim([0 length(w_hat)+1])
flabel = cellfun(@(f)functions(f).function,[features{:}],'Uni',0);
flabel = cellfun(@(s) s(strfind(s,')')+1:end),flabel,'Uni',0);
set(gca,'XTickLabel',flabel)
grid on

%%% p-values
figure(fignum)
fignum=fignum+1;
pvals = arrayfunvec(res,@(v)outn(@swtest,v,2),1);
pvals_0 = arrayfunvec(res_0,@(v)outn(@swtest,v,2),1);
plot(pvals,'o-')
xlim([1 size(w_hat_its,2)])
title(['p-val({OLS})=',num2str(pvals_0(1)),...
', p-val({WENDy})=',num2str(pvals(end))],'fontsize',9)
xlabel('iter')
legend('p-val({\bf w^{(n)}})','location','best','box','off')
grid on

%%% data
figure(fignum)
fignum=fignum+1;
h1=plot(1:length(tobs),xsub,'k-','linewidth',2);
hold on;
h2=plot(1:length(tobs),xobs,'r.','markersize',8);
xlim([1 length(tobs)])
if toggle_ddd
    h3=plot(1:length(t_learned),x_learned,'--g','linewidth',2);
    ylim([min(xobs(:)) max(xobs(:))+0.5*range(xobs(:))])
    try
        err_dd = norm(xsub(:) - x_learned(:))/norm(xsub(:));
    catch
        err_dd = Inf;
    end
    legend([h1(1);h2(1);h3(1)],{'u^*','{\bfU}',['{\bfU_{dd}}: ',num2str(100*err_dd),'% rel. err']},'location','best','box','off')
else
    legend([h1(1);h2(1)],{'u^*','{\bfU}'},'location','best','box','off')
end
title(['data: (K,d,M)=(',num2str(K),',',num2str(nstates),',',num2str(M),')'],'fontsize',9)
hold off
xlabel('timeindex')
grid on

%%% residual
figure(fignum)
fignum=fignum+1;
plot(norm(RT)*res(:,end),'LineWidth',2); 
xlim([1 length(b_0)])
legend({'{\bf w}_{WENDy}'},'box','off')
title(['\bf C^{-1/2}r(U,w). p-val ',num2str(pvals(end))],'fontsize',9)
xlabel('row num (k)')
grid on

%%% FFT data, FFT test function
figure(fignum)
fignum=fignum+1;
h1=semilogy(xfft(1:floor(end/2),:),'r');hold on; 
h2=semilogy(phifft(:,1:floor(end/2))'); hold off
ylim([min(min(xfft))/10 1])
xlim([1 length(xfft(1:floor(end/2),:))])
set(gca,'Ytick',10.^(floor(log10(min(min(xfft))/10)):0))
legend([h1(1);h2(1)],{'F({\bfU})','F(\Phi(1,:))'},'fontsize',10,'box','off')
title('Fourier content of data','fontsize',9)
xlabel('wavenumber')
grid on

%%% OLS residual: WENDy final to true
figure(fignum)
fignum=fignum+1;   
plot(res_0(:,1),'LineWidth',2)
hold off;
xlim([1 length(b_0)])
legend('{\bf w}_{WENDy}','box','off')
title(['\bf r(U,w). p-val ',num2str(outn(@swtest,res_0(:,1),2))],'fontsize',9)
xlabel('row num (k)')
if max(abs(res_0(:,1)))
    ylim(max(abs(res_0(:,1)))*[-1.5 1.5])
end
grid on




