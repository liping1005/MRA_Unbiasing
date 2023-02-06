figure
axis square
imagesc('XData',w,'YData',w,'CData',est_bispectrum_real)
title('Recovered','fontsize',16,'Interpreter','Latex')
colorbar()

figure
axis square
imagesc('XData',w,'YData',w,'CData',real(data))
title('Data Term','fontsize',16,'Interpreter','Latex')
colorbar()

figure
axis square
imagesc('XData',w,'YData',w,'CData',real(TrueBS))
title('Ground Truth','fontsize',16,'Interpreter','Latex')
colorbar()

figure
axis square
imagesc('XData',w,'YData',w,'CData',real(unbiased_BS))
title('unbiased BS','fontsize',16,'Interpreter','Latex')
colorbar()

figure
axis square
imagesc('XData',w,'YData',w,'CData',real(lowpass_BS))
title('low BS','fontsize',16,'Interpreter','Latex')
colorbar()