figure
subplot(1,3,1)
hold on
grid on
axis square
for i=1:min(M,10)
    plot(w, PowerSpectrum(i,:))
end
xlim([min(w) max(w)])
xlabel('$\omega$','fontsize',14,'Interpreter','Latex');
title('Target Power Spectrum and Power Spectrum for $M=10$ Observed Signals','fontsize',16,'Interpreter','Latex')
plot(w, UndilatedPowerSpectrum,'-.','Color',[0.4940, 0.1840, 0.5560],'LineWidth',2)

pos = get(gcf,'position');
set(gcf,'position',[6*pos(1:3) 2*pos(4)])

subplot(1,3,2)
plot(w,MeanPowerSpectrum-2*N*noise_sigma^2,'Linewidth',2,'Color',[0, 0.4470, 0.7410])
hold on
grid on
axis square
plot(w,TargetPowerSpectrum,'Linewidth',2,'Color',[0.4660 0.6740 0.1880])
xlim([min(w) max(w)])
xlabel('$\omega$','fontsize',14,'Interpreter','Latex');
title('Mean PS and Smoothed Mean PS','fontsize',16,'Interpreter','Latex')
legend({'$\widetilde{g}_\eta+\widetilde{g}_\sigma$','$(\widetilde{g}_\eta+\widetilde{g}_\sigma)\ast \phi_L$'},'FontSize',18,'Interpreter','latex')

subplot(1,3,3)
plot(w,UnbiasedPS,'Linewidth',2,'Color',[0.8500, 0.3250, 0.0980])
hold on
grid on
axis square
xlim([min(w),max(w)])
plot(w,MeanPowerSpectrum-2*N*noise_sigma^2,'Linewidth',2,'Color',[0, 0.4470, 0.7410])
plot(w,UndilatedPowerSpectrum,'-.','LineWidth',2,'Color',[0.4940, 0.1840, 0.5560])
xlabel('$\omega$','Fontsize',18,'Interpreter','latex')
title('PS Estimator from Inversion Unbiasing, Mean PS, and Target PS','fontsize',16,'Interpreter','Latex')
legend({'$\widetilde{Pf}$','$\widetilde{g}_\eta+\widetilde{g}_\sigma$','$Pf$'},'FontSize',18,'Interpreter','latex')
legend('Location','northeast')