% define spatial and freq resolution for bispectrum approx:

N = 2^4;
l = 2;
t = -(N):(1/2^l):(N)-1/2^l;
w = -pi*(2^l):(pi/N):pi*(2^l)-(pi/N);

%%
f = fftshift(exp(-t.^2/2)/sqrt(2*pi));
f2 = fftshift(exp(-(t-2).^2/2)/sqrt(2*pi));

f_hat = fft(f)*(1/2^l);
PS_f = f_hat.*conj(f_hat);

f2_hat = fft(f2)*(1/2^l);
PS_f2 = f2_hat.*conj(f2_hat);

Bf = compute_bispectrum3(f_hat);
Bf2 = compute_bispectrum3(f2_hat);

%Move blocks around to recenter frequencies, if needed:
Bf = center_BS(Bf);
Bf2 = center_BS(Bf2);

norm(real(Bf-Bf2))

norm(imag(Bf-Bf2))

norm(abs(Bf)-abs(Bf2))

norm(PS_f-PS_f2)

figure
imagesc('XData',w,'YData',w,'CData',real(Bf))
colorbar
title('original, real(Bf)')
set(gca,'YDir','normal')

figure
imagesc('XData',w,'YData',w,'CData',real(Bf2))
colorbar
title('shifted, real(Bf)')
set(gca,'YDir','normal')
% 
% figure
% imagesc('XData',w,'YData',w,'CData',imag(Bf))
% colorbar
% title('original, imag(Bf)')
% set(gca,'YDir','normal')
% 
% figure
% imagesc('XData',w,'YData',w,'CData',imag(Bf2))
% colorbar
% title('shifted, imag(Bf)')
% set(gca,'YDir','normal')


