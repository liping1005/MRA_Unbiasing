userpath('/mnt/ufs18/home-109/yinlipi1/bispectrum')
addpath(genpath('../../BS_inversion'))
addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../NumericalExperiments'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=4;
M=2^20;
t1=-(N/2):(1/2^l):(N/2)-1/2^l;
t = -(N):(1/2^l):(N)-1/2^l;
w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
spacing = w(2) - w(1);

Signal = 1; 
if Signal==1
    f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==2
    f1 = @(x)(10.6857).* exp(-5.*x.^2).*cos(12.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';    
elseif Signal==3
    f1 =@(x)1/sqrt(0.0078).*(sinc(4.*x));
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal == 4
    f1 = @(x) (8.1930* cos(6*x).*step_function(x, -.5, .5));
    RandomDilationOpts.SynthesisDomain = 'Space';
else 
    disp('Error')
end

f = [zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)];
Undilatedf_hat = ifftshift(fft(fftshift(f)));
UndilatedPowerSpectrum = abs(Undilatedf_hat).^2;
UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));

disp(sum(abs(UndilatedBispectrum).^2, 'all'))

