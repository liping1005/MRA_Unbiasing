% Playing around to see if we can construct medium freq example with decent
% norm on bispectrum

N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=4;

addpath(genpath('../Utils'))
addpath(genpath('../BS_Calculation'))



%f1 = @(x)exp(-20*x.^2);

%f1 = @(x)(step_function(x,-1,-.75)+step_function(x,-.25,.25)+step_function(x,.75,1))
f1 = @(x)(step_function(x,-.25,.25))

%f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);
%RandomDilationOpts.SynthesisDomain = 'Space';

%f1 = @(x)(10.6857)*exp(-20*x.^2).*cos(16.*x);
%RandomDilationOpts.SynthesisDomain = 'Space';

t1=-(N/2):(1/2^l):(N/2)-1/2^l;
t = -(N):(1/2^l):(N)-1/2^l;

% Pad with zeros:
f = [zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)];

% Define frequencies, switch to bispectrum here
w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
spacing = w(2) - w(1);
%get fft
Undilatedf_hat = fft(fftshift(f)).*(1/2^l);
UndilatedPowerSpectrum = abs(Undilatedf_hat).^2;
UndilatedBispectrum = CenterBS(ComputeBispectrum(Undilatedf_hat));

figure
axis square
imagesc('XData',w,'YData',w,'CData',real(UndilatedBispectrum))
title('Bispectrum','fontsize',16,'Interpreter','Latex')
colorbar()

% figure
% axis square
% imagesc('XData',w,'YData',w,'CData',imag(UndilatedBispectrum))
% title('Bispectrum','fontsize',16,'Interpreter','Latex')
% colorbar()

signal_norm = sqrt(sum(f.^2)*(t(2)-t(1)))

ft_norm = sqrt(sum(UndilatedPowerSpectrum)*(w(2)-w(1)))

ft_norm_sq = sum(UndilatedPowerSpectrum)*(w(2)-w(1))

BS_norm = sqrt(sum(sum(abs(UndilatedBispectrum).^2))*(w(2)-w(1))^2)

ratio = BS_norm/ft_norm_sq

% figure
% plot(t,f)