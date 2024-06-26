addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

%% Select example: 

Signal = input(['\nWhich example do you want to run?',...
       '\n ',...
       '\n1 Low Low Frequency Gabor',...
       '\n ',...
       '\n2 Low Frequency Gabor',...
       '\n ',...
       '\n3 Medium Frequency Gabor',...
       '\n ',...
       '\n4 Step Function',...
       '\n ',...
       '\n5 Squared Sinc',...
       '\n ',...
       '\n6 Triangle',...
       '\n ',...
       '\n7 Zero Function',...
       '\n ',...
       '\nInput example number without brackets or parentheses: ']);

if Signal==1
    f1 = @(x)(9.759)*exp(-5*x.^2) .* cos(4.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';  
elseif Signal==2
    f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==3 % consider making 16 the highest?
    f1 = @(x)(10.6857)*exp(-5*x.^2).*cos(12.*x);
    RandomDilationOpts.SynthesisDomain = 'Space'; 
elseif Signal==4 %step function is probably too hard
                 %bispectrum too small if space support large
                 %too tall and narrow if its space support small
    f1 = @(x)(5.659*step_function(x,-0.5,0.5));
    RandomDilationOpts.SynthesisDomain = 'Space';   
elseif Signal==5
    f1 = @(x)19.6116/2*(sinc(4.*x)).^2;
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==6
    f1 = @(x) 3.472*triangle(x,-1.0, 1.0, 2);
     RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==7
    f1 = @(x)zeros(size(x));
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==8
    f1 =@(x)1/sqrt(0.0078).*(sinc(4.*x));
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==9
    f1 = @(x)1/sqrt(0.0208).*(sinc(x)).^2;
else 
    disp('Error')
end


% Define parameters and signal (signal defined on [-N/2, N/2), noise on [-N,N) with spacing 1/2^l)
N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=3;
M=500000; % number of times we sample the noisy signal
GaussianConstant = 5.0;
PSWidthConstant = 5.0;
RandomDilationOpts.SynthesisDomain = 'Space';
RandomDilationOpts.Normalization = 'Linf'; %Options: L1 or Linf normalized dilations
RandomDilationOpts.Translate = 'True'; 
RandomDilationOpts.MagnitudeMaxTau = 0.25; %Needed for both Uniform and TruncatedGaussian (default: 0.2)
RandomDilationOpts.Distribution='Uniform';
RandomDilationOpts.UnbiasingMethod = 'Uniform';
RandomDilationOpts.PSUniform = 'yes'; %for UnbiasingMethod = 'Uniform'; options: yes or no
RandomDilationOpts.SmoothPS = 'yes';
RandomDilationOpts.SmoothDerivPS = 'no';
RandomDilationOpts.SmoothPSCorrectionTerm = 'no'; %Only compute for Oracle! Hasn't been implemented for Empirical
RandomDilationOpts.InterpolationMethod = 'spline';
RandomDilationOpts.MomentCalc='Empirical';
RandomDilationOpts.DilationCalc=RandomDilationOpts.MomentCalc; 
true_noise_sigma = 1.0; %Optional: add additive Gaussian noise
OptimizationOpts.Method = 'Unconstrained';
OptimizationOpts.Initialization = 'MeanPS_NoDilUnbias'; %options: MeanPS_NoDilUnbias, MeanPS_Order2Unbias, MeanPS_Order4Unbias
OptimizationOpts.tol = 1e-7;
OptimizationOpts.Uniformtol = 1e-7; % when UnbiasingMethod = 'Uniform', tolerance for recovering g from h via optimization
PlotFigs = 'no'; %options: 'yes' or 'no'

ProcessSignalsV2
eta_true = sqrt(var(Tau));
estimated_noise_sigma = sqrt(mean([PowerSpectrumAvg(1:N) PowerSpectrumAvg(end-N+1:end)])/(2*N));
if strcmp(RandomDilationOpts.MomentCalc,'Oracle')
    eta = eta_true;
    if strcmp(RandomDilationOpts.SmoothPS, 'yes')
        MakeSmoothingMatrixPS
        TargetPowerSpectrum = (SmoothingMatrix*(PowerSpectrumAvg - 2*N*estimated_noise_sigma^2)')';
    elseif strcmp(RandomDilationOpts.SmoothPS, 'no')
        TargetPowerSpectrum = PowerSpectrumAvg - 2*N*estimated_noise_sigma^2;
    end
    if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
        Unbias_Uniform
    elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
        Unbias_GeneralDerivative
    end
elseif strcmp(RandomDilationOpts.MomentCalc,'Empirical')
    if strcmp(RandomDilationOpts.SmoothPS, 'yes')
        MakeSmoothingMatrixPS
        TargetPowerSpectrum = (SmoothingMatrix*(PowerSpectrumAvg - 2*N*estimated_noise_sigma^2)')';
    elseif strcmp(RandomDilationOpts.SmoothPS, 'no')
        TargetPowerSpectrum = PowerSpectrumAvg - 2*N*estimated_noise_sigma^2;
    end
    if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
        Unbias_Uniform
    elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
        Unbias_GeneralDerivative
    end
    eta = eta_PS;
    disp('True Eta:')
    disp(eta_true) 
    disp('Estimated Eta:')
    disp(eta)  
end

figure
plot(w,UnbiasedPS)
hold on
plot(w,UndilatedPowerSpectrum)
plot(w,TargetPowerSpectrum)
legend({'Unbiased','True','MeanPS'},'fontsize',14)