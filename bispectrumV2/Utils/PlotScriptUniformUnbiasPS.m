plot_idx = 1:1:length(Mvalues);
errorbar(log2(Mvalues(plot_idx)), log2(PSerrorVec_NoDilUnbias(plot_idx)),SDLogPSerrorVec_NoDilUnbias(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
xlim([min(log2(Mvalues(plot_idx))) max(log2(Mvalues(plot_idx)))])

axis square
grid on
hold on
errorbar(log2(Mvalues(plot_idx)), log2(PSerrorVec(plot_idx)), SDLogPSerrorVec(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
legend({'PS (No Dilation UB)', 'PS (Inversion UB)'},'Fontsize',12,'Location','Southwest')
xlabel('$\log_2(M)$','Fontsize',18,'Interpreter','latex')
ylabel('$\log_2$(Error)','Fontsize',18,'Interpreter','latex')