f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);

N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=3;
t1=-(N/2):(1/2^l):(N/2)-1/2^l;
t = -(N):(1/2^l):(N)-1/2^l;
w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
spacing = w(2) - w(1);


f = [zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)];
f_noisy = f + randn(size(f))/25.0;
f_hat_noisy = ifftshift(fft(fftshift(f_noisy)));
ps_noisy = abs(f_hat_noisy);
bs_noisy = CenterBS(ComputeBispectrum(fftshift(f_hat_noisy)));


invert_bispectrum(f_hat_noisy, ps_noisy, bs_noisy,f, w, t, l)
