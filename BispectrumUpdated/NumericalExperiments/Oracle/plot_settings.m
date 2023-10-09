N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=5;
Mvalues = 2.^(4:1:20);
NumberSimulationsPerValue = 5;

GaussianConstant = 5;
PSWidthConstant = 10;
t1=-(N/2):(1/2^l):(N/2)-1/2^l;
w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
spacing = w(2) - w(1);

RandomDilationOpts.Normalization = 'Linf'; %Options: L1 or Linf normalized dilations
RandomDilationOpts.Translate = 'True'; 
RandomDilationOpts.MagnitudeMaxTau = 0.5; %Needed for both Uniform and TruncatedGaussian (default: 0.2)
RandomDilationOpts.Distribution='Uniform';
RandomDilationOpts.UnbiasingMethod = 'Uniform';
RandomDilationOpts.PSUniform = 'yes'; %for UnbiasingMethod = 'Uniform'; options: yes or no
RandomDilationOpts.SmoothPS = 'yes';
RandomDilationOpts.SmoothDerivPS = 'no';
RandomDilationOpts.SmoothPSCorrectionTerm = 'no'; %Only compute for Oracle! Hasn't been implemented for Empirical
RandomDilationOpts.InterpolationMethod = 'linear';
RandomDilationOpts.DilationCalc='Empirical';
true_noise_sigma = 1.0; %Optional: add additive Gaussian noise
%sqrt(2)
OptimizationOpts.Method = 'Unconstrained';
OptimizationOpts.Initialization = 'MeanPS_NoDilUnbias'; %options: MeanPS_NoDilUnbias, MeanPS_Order2Unbias, MeanPS_Order4Unbias
OptimizationOpts.tol = 1e-7;
OptimizationOpts.Uniformtol = 1e-7; % when UnbiasingMethod = 'Uniform', tolerance for recovering g from h via optimization
PlotFigs = 'no'; %options: 'yes' or 'no'

