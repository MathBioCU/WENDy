clf
semilogy(1:size(w_hat_its,2)-1,vecnorm(diff(w_hat_its,[],2))/norm(true_vec),'b*-',...
    0:size(w_hat_its,2)-1,vecnorm(w_hat_its-w_hat)/norm(true_vec),'ro-',...
    0:size(w_hat_its,2)-1,vecnorm(w_hat_its-true_vec)/norm(true_vec),'gd-','linewidth',2,'markersize',8)
legend({'||w^{(n+1)}-w^{(n)}||/||w^*||','||w^{(n+1)}-w^{\^} ||/||w^*||','||w^{(n+1)}-w^*||/||w^*||'})
grid on
xlabel('iteration')
ylabel('log_{10}(error)')
set(gca,'fontsize',16)
xlim([1,size(w_hat_its,2)-1])

axes('position',[0.2 0.17 0.32 0.3])
% h1=plot(1:length(tobs),xsub,'k-','linewidth',2);
% hold on;
% h2=plot(1:length(tobs),xobs,'r.','markersize',8);
% xlim([1 length(tobs)])
% if toggle_ddd
%     h3=plot(1:length(t_learned),x_learned,'--g','linewidth',2);
%     ylim([min(xobs(:)) max(xobs(:))+0.5*range(xobs(:))])
%     try
%         err_dd = norm(xsub(:) - x_learned(:))/norm(xsub(:));
%     catch
%         err_dd = Inf;
%     end
%     legend([h1(1);h2(1);h3(1)],{'u^*','{\bfU}',['{\bfU_{dd}}: ',num2str(100*err_dd),'% rel. err']},'location','best','box','off')
% else
%     legend([h1(1);h2(1)],{'u^*','{\bfU}'},'location','best','box','off')
% end
% title(['data: (K,d,M)=(',num2str(K),',',num2str(nstates),',',num2str(M),')'],'fontsize',9)
% hold off
% xlabel('timeindex')
% grid on

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
xlabel('parameter index')
grid on

% saveas(gcf,['~/Dropbox/Boulder/research/writing/research_log/figures/wendy_conv_HR_02nz_1024M.png'])
