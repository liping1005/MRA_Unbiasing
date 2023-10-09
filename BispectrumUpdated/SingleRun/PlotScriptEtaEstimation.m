Signal = 15;
l = 4;
eta_true = 12^(-0.5);
myStr = sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/Empirical/results_f%d_l%d_learn_eta.mat', Signal, l);
load(myStr)
gcf = figure
Eta_Error_AllSims = abs(Eta_PS_AllSims - eta_true)/eta_true
MeanEtaError = mean(Eta_Error_AllSims,2)
SDEtaError = std(Eta_Error_AllSims')'
errorbar(log2(Mvalues), log2(MeanEtaError),log2(1-SDEtaError./(MeanEtaError.*sqrt(NumberSimulationsPerValue))),log2(1+SDEtaError./(MeanEtaError.*sqrt(NumberSimulationsPerValue))),'-','Linewidth',2)
xlim([min(log2(Mvalues(plot_idx))) max(log2(Mvalues(plot_idx)))])
ylim([-8, 0])
xlabel('$\log_2(M)$','Fontsize',18,'Interpreter','latex')
ylabel('$\log_2$($\eta$-Error)','Fontsize',18,'Interpreter','latex')
axis square
grid on

save_name = sprintf('f%d_eta_error.pdf', Signal);
exportgraphics(gcf, save_name)