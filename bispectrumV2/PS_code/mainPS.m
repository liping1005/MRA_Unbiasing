addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

%% Select example: 

Signal = input(['\nWhich example do you want to run?',...
       '\n ',...
       '\n1 Low Frequency Gabor',...
       '\n ',...
       '\n2 Medium Frequency Gabor',...
       '\n ',...
       '\n3 High Frequency Gabor',...
       '\n ',...
       '\n4 Sinc in Frequency',...
       '\n ',...
       '\n5 High Frequency Chirp',...
       '\n ',...
       '\n6 Step in Frequency',...
       '\n ',...
       '\n7 Zigzag in Frequency',...
       '\n ',...
       '\nInput example number without brackets or parentheses: ']);

if Signal==1
    
    f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
    
elseif Signal==2
    
    f1 = @(x)(10.6857)*exp(-5*x.^2).*cos(16.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
    
elseif Signal==3
    
    f1 = @(x)(10.6857)*exp(-5*x.^2).*cos(32.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
    
elseif Signal==4
    
    f1 = @(x)(4.45458)*(sinc(.2.*(x-32))+sinc(.2.*(-x-32)));
    RandomDilationOpts.SynthesisDomain = 'Frequency';
    
elseif Signal==5
    
    f1 = @(x)(3.19584)*exp(-.04*(x).^2).*cos(30*(x)+1.5*x.^2);
    RandomDilationOpts.SynthesisDomain = 'Space';   
    
elseif Signal==6
    
    f1 = @(x)(4.09331)*(step_function(x,-38,-32)+step_function(x,32,38));
    RandomDilationOpts.SynthesisDomain = 'Frequency';
    
elseif Signal==7
    
    f1 = @(x)(2.58883)*sqrt(zigzag((x+40)/5)+zigzag((x-40)/5));
    RandomDilationOpts.SynthesisDomain = 'Frequency';
else 
    disp('Error')
end


% Define parameters and signal (signal defined on [-N/2, N/2), noise on [-N,N) with spacing 1/2^l)
N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=5;
M=5000; % number of times we sample the noisy signal
RandomDilationOpts.SynthesisDomain = 'Space';
RandomDilationOpts.Normalization = 'Linf'; %Options: L1 or Linf normalized dilations
RandomDilationOpts.Translate = 'True'; 
RandomDilationOpts.MagnitudeMaxTau = 0.5; %Needed for both Uniform and TruncatedGaussian (default: 0.2)
RandomDilationOpts.Distribution='Uniform';
RandomDilationOpts.UnbiasingMethod = 'Uniform';
RandomDilationOpts.PSUniform = 'yes'; %for UnbiasingMethod = 'Uniform'; options: yes or no
RandomDilationOpts.SmoothPS = 'yes';
RandomDilationOpts.SmoothDerivPS = 'no';
RandomDilationOpts.SmoothPSCorrectionTerm = 'no'; %Only compute for Oracle! Hasn't been implemented for Empirical
RandomDilationOpts.InterpolationMethod = 'spline';
RandomDilationOpts.MomentCalc='Empirical';
%RandomDilationOpts.MomentCalc='Empirical';
true_noise_sigma = sqrt(2); %Optional: add additive Gaussian noise
OptimizationOpts.Method = 'Unconstrained';
OptimizationOpts.Initialization = 'MeanPS_NoDilUnbias'; %options: MeanPS_NoDilUnbias, MeanPS_Order2Unbias, MeanPS_Order4Unbias
OptimizationOpts.tol = 1e-7;
OptimizationOpts.Uniformtol = 1e-7; % when UnbiasingMethod = 'Uniform', tolerance for recovering g from h via optimization
PlotFigs = 'yes'; %options: 'yes' or 'no'

ProcessSignalsV2
eta_true = sqrt(var(Tau));
estimated_noise_sigma = sqrt(mean([PowerSpectrumAvg(1:N) PowerSpectrumAvg(end-N+1:end)])/(2*N));
if strcmp(RandomDilationOpts.MomentCalc,'Oracle')
    eta = eta_true;
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