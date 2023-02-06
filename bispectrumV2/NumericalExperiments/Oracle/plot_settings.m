N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=5;
%Mvalues = ceil(2.^(17:.5:20)); %increasing M
Mvalues = 2.^(4:2:16);
%Mvalues = 2.^(8:2:12);
NumberSimulationsPerValue = 5;

RandomDilationOpts.Normalization = 'Linf'; %Options: L^1 or Linf normalized dilations
%RandomDilationOpts.Normalization = 'L1';
RandomDilationOpts.Translate = 'True'; %Only use 'True' when SynthesisDomain = 'Space'!
RandomDilationOpts.SynthesisDomain = 'Space'; 
%RandomDilationOpts.SynthesisDomain = 'Frequency'; 
RandomDilationOpts.MagnitudeMaxTau = 0.5; %Needed for both Uniform and TruncatedGaussian (default: 0.2)
RandomDilationOpts.Distribution='Uniform';
%RandomDilationOpts.Distribution='TruncatedGaussian'; RandomDilationOpts.SDTruncatedGaussian = 2^(-2);
%RandomDilationOpts.Distribution='NoDilation';
%RandomDilationOpts.UnbiasingMethod = 'GeneralDerivative'; %options: 'GeneralDerivative' or 'Uniform'
RandomDilationOpts.UnbiasingMethod = 'Uniform';
RandomDilationOpts.PSUnbiasingOrder = 0; %for UnbiasingMethod = 'GeneralDerivative'; options: 0,2,4 %Highest even moment to use in unbiasing procedure
RandomDilationOpts.PSUniform = 'yes'; %for UnbiasingMethod = 'Uniform'; options: yes or no
RandomDilationOpts.SmoothPS = 'yes';
RandomDilationOpts.SmoothDerivPS = 'no';
RandomDilationOpts.SmoothPSCorrectionTerm = 'no';
RandomDilationOpts.InterpolationMethod = 'spline';
%RandomDilationOpts.DilationCalc='Oracle';
RandomDilationOpts.DilationCalc='Oracle';
true_noise_sigma = sqrt(2); %Optional: add additive Gaussian noise
OptimizationOpts.Method = 'Unconstrained';
OptimizationOpts.Initialization = 'MeanPS_NoDilUnbias'; %options: MeanPS_NoDilUnbias, MeanPS_Order2Unbias, MeanPS_Order4Unbias
OptimizationOpts.tol = 1e-7;
OptimizationOpts.Uniformtol = 1e-7; % when UnbiasingMethod = 'Uniform', tolerance for recovering g from h via optimization
PlotFigs = 'yes'; %options: 'yes' or 'no'

