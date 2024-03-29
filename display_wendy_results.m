%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

%% print results

disp(['-----------------'])
disp([' '])
disp('parameters: true, WENDy, NLS')
disp([true_vec w_hat param_nls])
disp(['err WENDy:',num2str(err_wendy)])
disp(['time WENDy: ',num2str(total_time), '(sec)'])
disp(['Num its: ', num2str(size(w_hat_its,2))])
disp(['K,d,M:',num2str(size(V_cell{1},1)),',',num2str(length(x0)),',',num2str(length(tobs))])
disp(['test function radii: (',num2str(unique(mt)'),')'])
disp([' '])
disp(['time NLS:',num2str(total_time_nls),'(sec)'])
disp(['err NLS:',num2str(err_NLS)])
disp([' '])
disp(['-----------------'])

%% compute quantities of interest

%%% compute (G*,b*)
x_cell = mat2cell(xsub,M,ones(1,nstates));
Theta_cell_true = cellfun(@(x) cell2mat(cellfun(@(y) y(x_cell{:}), x, 'uni',0)), features,'uni',0);
G_0_true = blkdiag(V_cell{:})*blkdiag(Theta_cell_true{:});
b_0_true = -blkdiag(Vp_cell{:})*xsub(:);    
K = size(b_0,1);

%%% Phi'*ep
response_error = b_0_true - b_0;

%%% Phi(Theta(U)-Theta(U-ep))
F_obs_error = G_0*w_hat - G_0_true*w_hat;

%%% Phi*u' + Phi'*u
int_error = (RT/norm(RT)) \ (G_0_true*true_vec - b_0_true);

%%% G*(w-w*)
w_error_response = (RT/norm(RT)) \ G_0_true*(w_hat - true_vec);

%%% Phi*nabla(F)*eps
J = get_jac(Jac_mat,w_hat,features);
Jac_F_obs_error = blkdiag(V_cell{:})*reshape(squeeze(pagemtimes(J,permute(noise,[2 3 1])))',[],1);
lin_approx = (RT/norm(RT)) \ (Jac_F_obs_error + response_error);

%%% nonlinear portion of residual
nonlin_res = (RT/norm(RT)) \ (F_obs_error - Jac_F_obs_error);

%%% full residual (same as res(:,end) up to scaling)
Res_full = w_error_response + int_error + lin_approx + nonlin_res;

%%% compute data-driven dynamics
if toggle_ddd
    ind=inf; % choose which weights to generate from
    tol_ode = 1e-12;
    w_jac_cell = vec2cell(w_hat_its(:,min(ind,end)),cellfun(@(x)length(x),features));
    options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
    [t_learned,x_learned]=ode45(@(t,x)rhs_p(x,w_jac_cell),tobs,x0,options_ode_sim);
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

err_all = vecnorm([true_vec w_hat param_nls]-true_vec)/norm(true_vec);

if toggle_plot
    figure(11)
    
    %%% wendy iterates
    subplot(3,3,1)
    semilogy(1:length(errs),errs,'bo-')
    legend({'err({\bf w^{(n)}})'},'location','best','box','off')
    ylabel('||\bf w^{(n)}-w^*||_2/||w^*||_2')
    title(['err({OLS})=',num2str(100*errs(1)),'%, ',...
        'err({WENDy})=',num2str(100*errs(end)),'%'],'fontsize',9)
    ylims = [10^(log10(min(errs)))*0.9 10^(log10(max(errs)))*1.2];
    ylim(ylims)
    xlim([1 length(errs)])
    set(gca,'Ytick',exp(linspace(log(ylims(1)),log(ylims(2)),8)))
    xlabel('iter')
    grid on

    %%% confidence intervals
%     subplot(3,3,2)
%     h=fill(xflip,confbounds,'red','FaceColor',[1 0.8 0.8],'edgecolor','none'); hold on;
%     % h1=plot(1:length(w_hat),(w_hat-true_vec)./abs(true_vec),'ro',1:length(w_hat),w_hat*0,'b--');
%     h1=plot(1:length(w_hat),w_hat,'ro',1:length(w_hat),true_vec,'bx');
%     hold off
%     legend(h1(1:2),{'{\bf w}_{WENDy}','\bf w^*'},'box','off','location','best')
%     title([num2str((1-c)*100),'% confidence bounds'],'fontsize',9)
%     set(gca,'Xtick',1:length(w_hat))
%     xlim([0 length(w_hat)+1])
%     grid on
% 
    subplot(3,3,2)
h1=plot(1:length(w_hat),w_hat,'ro',1:length(w_hat),true_vec,'bx');hold on
for j=1:length(true_vec)
    rectangle('position',[j-0.25 w_hat(j)-conf_int(j) 0.5 2*conf_int(j)]);
    line([j-0.25 j+0.25],[w_hat(j) w_hat(j)])
end
hold off
legend(h1(1:2),{'{\bf w}_{WENDy}','\bf w^*'},'box','off','location','best')
title([num2str((1-c)*100),'% confidence bounds'],'fontsize',9)
set(gca,'Xtick',1:length(w_hat))
xlim([0 length(w_hat)+1])
grid on

    % ylim([-1 1]*2*range((w_hat-true_vec)./abs(true_vec)))
    
    %%% p-values
    subplot(3,3,3)
    pvals = arrayfunvec(res,@(v)outn(@swtest,v,2),1);
    pvals_0 = arrayfunvec(res_0,@(v)outn(@swtest,v,2),1);
    plot(pvals,'o-')
    xlim([1 length(errs)])
    title(['p-val({OLS})=',num2str(pvals_0(1)),...
        ', p-val({WENDy})=',num2str(pvals(end))],'fontsize',9)
    xlabel('iter')
    legend('p-val({\bf w^{(n)}})','location','best','box','off')
    grid on

    %%% data
    subplot(3,3,4)
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
    
    %%% GLS residual: WENDy final to true
    subplot(3,3,5)
    plot(norm(RT)*res(:,end),'LineWidth',2); 
    hold on; 
    plot(norm(RT)*res_true(:,end),'r--','LineWidth',2); 
    hold off;
    xlim([1 length(b_0)])
    legend({'{\bf w}_{WENDy}','{\bf w}^*'},'box','off')
    title(['\bf C^{-1/2}r(U,w). p-val ',num2str(pvals(end))],'fontsize',9)
    xlabel('row num (k)')
    if max(abs(Res_full))
        ylim(max(abs(Res_full))*[-1.5 1.5])
    end
    grid on

    
    %%% residual components
    subplot(3,3,6)
    plot(1:K,Res_full,'k-','LineWidth',2)
    hold on
    plot(1:K,lin_approx,'r--','LineWidth',2);
    xlim([1 length(b_0)])
    hold off
%     title('residual components','fontsize',9)
    title('\bf C^{-1/2}r(U,w) vs. C^{-1/2}L_w\epsilon','fontsize',9)
    xlabel('row num (k)')
    legend({'\bfr(U,w)','\bfL_w\epsilon'},'location','best','box','off')
    if max(abs(Res_full))
        ylim(max(abs(Res_full))*[-1.5 1.5])
    end
    grid on


    %%% FFT data, FFT test function
    subplot(3,3,7)
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
    subplot(3,3,8)    
    plot(res_0(:,1),'LineWidth',2);hold on; 
    plot(res_0_true(:,1),'r--','LineWidth',2)
    hold off;
    xlim([1 length(b_0)])
    legend({'{\bf w}_{WENDy}','{\bf w}^*'},'box','off')
    title(['\bf r(U,w). p-val ',num2str(outn(@swtest,res_0(:,1),2))],'fontsize',9)
    xlabel('row num (k)')
    if max(abs(res_0(:,1)))
        ylim(max(abs(res_0(:,1)))*[-1.5 1.5])
    end
    grid on

    
    %%% residual and error of linear approx
    subplot(3,3,9)
    plot(1:length(Res_full),[w_error_response int_error nonlin_res],'-.','LineWidth',2.5);
    xlim([1 length(b_0)])
    if max(abs(Res_full))
        ylim(max(abs(Res_full))*[-1.5 1.5])
    end
    legend({'{\bfr}_0','{\bfe}_{int}','\bfh'},'location','best','box','off')
    title('other residual components','fontsize',9)
    xlabel('row num (k)')
    grid on


end
