real_data = real(data);
imag_data = imag(data);
real_TrueBS = real(UndilatedBispectrum);
imag_TrueBS = imag(UndilatedBispectrum);
real_lowpass = real(lowpass_BS);
imag_lowpass = imag(lowpass_BS);
%save('f4_real_C_5.mat', 'est_bispectrum_real', 'real_data', 'real_TrueBS','real_lowpass')
%save('f4_imag_C_5.mat', 'est_bispectrum_im', 'imag_data', 'imag_TrueBS', 'imag_lowpass')

figure
subplot(1,3,1)
hold on
grid on
axis square
imagesc('XData',w,'YData',w,'CData',real(Plotrecovered_bs))
title('Recovered','fontsize',16,'Interpreter','Latex')
colorbar()

subplot(1,3,2)
axis square
imagesc('XData',w,'YData',w,'CData',real(lowpass_BS))
title('BS After Centering','fontsize',16,'Interpreter','Latex')
colorbar()

subplot(1,3,3)
axis square
imagesc('XData',w,'YData',w,'CData',real(UndilatedBispectrum))
title('Ground Truth','fontsize',16,'Interpreter','Latex')
colorbar()

hold off

