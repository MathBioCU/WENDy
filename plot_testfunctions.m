%%% plot orthonormal test functions. Requires V_cell, Vp_cell, tobs
inds=1:6;
toggle_savefig = 1;
for i=1:length(inds)
    figure(i);clf
    set(gcf,'position',[1500+floor(i/4)*1000 1200-325*(mod(i-1,3)+1) 687 213])

    subplot(1,2,1)% plot tf
    plot(V_cell{1}(i,:)/mean(diff(tobs)),'linewidth',3)
    xlabel('timepoints')
    xlim([1 length(tobs)])
    ylim(max(max(Vp_cell{1}(inds,:)/mean(diff(tobs))))*[-1.5 1.5])
    grid on
    legend(strcat('\phi_{',num2str(inds(i)),'}'),'fontsize',12)
    set(gca,'fontsize',12)

    subplot(1,2,2)% plot tf deriv
    plot(Vp_cell{1}(i,:)/mean(diff(tobs)),'linewidth',3)    
    xlabel('timepoints')        
    grid on
    xlim([1 length(tobs)])
    ylim(max(max(Vp_cell{1}(inds,:)/mean(diff(tobs))))*[-1.5 1.5])
    legend(strcat('\phi^\prime_{',num2str(inds(i)),'}'),'fontsize',12)
    set(gca,'fontsize',12)

    if toggle_savefig 
        saveas(gcf,['~/Desktop/phi',strrep(num2str(inds(i)),' ',''),'.png'])
    end
end