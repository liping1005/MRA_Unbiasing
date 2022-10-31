figure
hold on
grid on
axis square
imagesc('XData',w,'YData',w,'CData',abs(est_bispectrum_real))
title('Real Estimate','fontsize',16,'Interpreter','Latex')
colorbar()

% figure
% hold on
% grid on
% axis square
% imagesc('XData',w,'YData',w,'CData',est_bispectrum_im)
% title('Imaginary Estimate','fontsize',16,'Interpreter','Latex')
% colorbar()

figure
hold on
grid on
axis square
imagesc('XData',w,'YData',w,'CData',real(BispectrumAvg))
title('Real Noisy','fontsize',16,'Interpreter','Latex')
colorbar()

% 
% figure
% hold on
% grid on
% axis square
% imagesc('XData',w,'YData',w,'CData',imag(BispectrumAvg))
% title('Imaginary Noisy','fontsize',16,'Interpreter','Latex')
% colorbar()

figure
hold on
grid on
axis square
imagesc('XData',w,'YData',w,'CData',real(TrueBS))
title('Real Actual','fontsize',16,'Interpreter','Latex')
colorbar()

% figure
% hold on
% grid on
% axis square
% imagesc('XData',w,'YData',w,'CData',imag(TrueBS))
% title('Imaginary Actual','fontsize',16,'Interpreter','Latex')
% colorbar()