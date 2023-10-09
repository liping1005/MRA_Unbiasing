Signal = 15;
l = 4;
myStr = sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/Empirical/results_f%d_l%d_sigma_half_learn_eta.mat', Signal, l);
load(myStr);

plot_idx = 1:1:length(Mvalues);
gcf = figure;
axis square
grid on
hold on
errorbar(log2(Mvalues(plot_idx)), log2(BSerrorVec_NoDilUnbias(plot_idx)),SDLogBSerrorVec_NoDilUnbias(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
xlim([min(log2(Mvalues(plot_idx))) max(log2(Mvalues(plot_idx)))])
[slope_no_ub_half, ~] = polyfit(log2(Mvalues),log2(BSerrorVec_NoDilUnbias),1);
errorbar(log2(Mvalues(plot_idx)), log2(BSerrorVec(plot_idx)), SDLogBSerrorVec(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
[slope_ub_half, ~] = polyfit(log2(Mvalues),log2(BSerrorVec), 1);
myStr = sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/Empirical/results_f%d_l%d_learn_eta.mat', Signal, l);
load(myStr)
errorbar(log2(Mvalues(plot_idx)), log2(BSerrorVec_NoDilUnbias(plot_idx)),SDLogBSerrorVec_NoDilUnbias(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
[slope_no_ub, ~] = polyfit(log2(Mvalues),log2(BSerrorVec_NoDilUnbias), 1);
errorbar(log2(Mvalues(plot_idx)), log2(BSerrorVec(plot_idx)), SDLogBSerrorVec(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
[slope_ub ,~]= polyfit(log2(Mvalues),log2(BSerrorVec), 1);
label_no_ub_half = sprintf('sigma = 0.5 (No UB)[Slope = %.2f]', slope_no_ub_half(1));
label_ub_half= sprintf('sigma = 0.5 (UB)[Slope = %.2f]', slope_ub_half(1));
label_no_ub= sprintf('sigma = 1.0 (No UB)[Slope = %.2f]', slope_no_ub(1));
label_ub= sprintf('sigma = 1.0 (UB)[Slope = %.2f]', slope_ub(1));
legend({label_no_ub_half, label_ub_half, label_no_ub, label_ub},'Fontsize',12,'Location','northeast')
xlabel('$\log_2(M)$','Fontsize',18,'Interpreter','latex')
ylabel('$\log_2$(Error)','Fontsize',18,'Interpreter','latex')
hold off

hold off
save_name = sprintf('f%d_learn_eta_bs_recovery.pdf', Signal);
exportgraphics(gcf, save_name)