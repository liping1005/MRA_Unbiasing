Signal = 2;
l = 4;
true_noise_sigma = 1.0;
if(true_noise_sigma == 0.5)
    myStr = sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/Empirical/results_f%d_l%d_oracle_sigma_half.mat', Signal, l);
elseif(true_noise_sigma == 1.0)
    myStr = sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/Empirical/results_f%d_l%d_learn_eta.mat', Signal, l);
end

load(myStr);

plot_idx = 1:1:length(Mvalues);
gcf = figure;
axis square
grid on
hold on
errorbar(log2(Mvalues(plot_idx)), log2(SignalerrorVec_NoDilUnbias_APS(plot_idx)),SDLogSignalerrorVec_NoDilUnbias_APS(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
xlim([min(log2(Mvalues(plot_idx))) max(log2(Mvalues(plot_idx)))])
[slope_no_ub_half ,~]= polyfit(log2(Mvalues),log2(SignalerrorVec_NoDilUnbias_APS), 1);
errorbar(log2(Mvalues(plot_idx)), log2(SignalerrorVec_APS(plot_idx)), SDLogSignalerrorVec_APS(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
[slope_ub_half ,~]= polyfit(log2(Mvalues),log2(SignalerrorVec_APS), 1);
errorbar(log2(Mvalues(plot_idx)), log2(SignalerrorVec_NoDilUnbias_FM(plot_idx)),SDLogSignalerrorVec_NoDilUnbias_FM(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
[slope_no_ub ,~]= polyfit(log2(Mvalues),log2(SignalerrorVec_NoDilUnbias_FM), 1);
errorbar(log2(Mvalues(plot_idx)), log2(SignalerrorVec_FM(plot_idx)), SDLogSignalerrorVec_FM(plot_idx)/sqrt(NumberSimulationsPerValue),'-','Linewidth',2)
[slope_ub ,~]= polyfit(log2(Mvalues),log2(SignalerrorVec_FM), 1);
label_no_ub_half = sprintf('APS (NO UB) [Slope = %.2f]', slope_no_ub_half(1));
label_ub_half= sprintf('APS (UB)[Slope = %.2f]', slope_ub_half(1));
label_no_ub= sprintf('FM (No UB)[Slope = %.2f]', slope_no_ub(1));
label_ub= sprintf('FM (UB)[Slope = %.2f]', slope_ub(1));
xlabel('$\log_2(M)$','Fontsize',18,'Interpreter','latex')
ylabel('$\log_2$(Error)','Fontsize',18,'Interpreter','latex')
legend({label_no_ub_half, label_ub_half, label_no_ub, label_ub},'Fontsize',12,'Location','northeast')
if(true_noise_sigma ==0.5)
    save_name = sprintf('f%d_inv_sigma_half.pdf', Signal);
elseif(true_noise_sigma == 1.0)
    save_name = sprintf('f%d_inv.pdf', Signal);
end
exportgraphics(gcf, save_name)