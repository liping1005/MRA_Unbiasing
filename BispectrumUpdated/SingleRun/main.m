addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

%% Set Parameters  

N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=3;
M= 50000; % number of times we sample the noisy signal
GaussianConstant = 1;
PSWidthConstant = 1;
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
RandomDilationOpts.DilationCalc='Oracle';
true_noise_sigma = 0.0; %Optional: add additive Gaussian noise
%sqrt(2)
OptimizationOpts.Method = 'Unconstrained';
OptimizationOpts.Initialization = 'MeanPS_NoDilUnbias'; %options: MeanPS_NoDilUnbias, MeanPS_Order2Unbias, MeanPS_Order4Unbias
OptimizationOpts.tol = 1e-7;
OptimizationOpts.Uniformtol = 1e-7; % when UnbiasingMethod = 'Uniform', tolerance for recovering g from h via optimization
PlotFigs = 'no'; %options: 'yes' or 'no'

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
    f1 = @(x)(9.759)*exp(-5*x.^2).*cos(4.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';  
elseif Signal==2
    f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==3 % consider making 16 the highest?
    f1 = @(x)(10.6857)*exp(-5*x.^2).*cos(16.*x);
    RandomDilationOpts.SynthesisDomain = 'Space'; 
elseif Signal==4
    f1 = @(x)(4.0 * step_function(x,-1.0,1.0));
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
else 
    disp('Error')
end


ProcessSignals

figure
plot(w,PlotUnbiasedPS)
hold on
plot(w,UndilatedPowerSpectrum)
plot(w,PlotTargetPowerSpectrum)
legend({'Unbiased','True','MeanPS'},'fontsize',14)
hold off

PlotResults

inversion