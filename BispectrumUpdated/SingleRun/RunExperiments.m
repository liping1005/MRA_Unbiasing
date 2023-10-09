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

GaussianConstant = 5;
PSWidthConstant = 5;
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
RandomDilationOpts.MomentCalc = RandomDilationOpts.DilationCalc;
true_noise_sigma = 0.5; %Optional: add additive Gaussian noise
%sqrt(2)
OptimizationOpts.Method = 'Constrained';
OptimizationOpts.Initialization = 'MeanPS_NoDilUnbias'; %options: MeanPS_NoDilUnbias, MeanPS_Order2Unbias, MeanPS_Order4Unbias
OptimizationOpts.tol = 1e-5;
OptimizationOpts.Uniformtol = 1e-5; % when UnbiasingMethod = 'Uniform', tolerance for recovering g from h via optimization
PlotFigs = 'no'; %options: 'yes' or 'no'

t1=-(N/2):(1/2^l):(N/2)-1/2^l;
t = -(N):(1/2^l):(N)-1/2^l;
w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
spacing = w(2) - w(1);

Signal = 2;
 
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

%load data
if(true_noise_sigma == 0.5)
    myStr =  sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/generatedSignals/f%d_signals_all_l%d_sigma_half.mat', Signal, l);
elseif(true_noise_sigma == 1.0)
    myStr =  sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/generatedSignals/f%d_signals_all_l%d.mat', Signal, l);   
end 
load(myStr)

FourierTransformAvg = reshape(FTAvgArray(1,3,9,:), [1,size(FTAvgArray, 4)]); 
PowerSpectrumAvg = reshape(PSAvgArray(1,3,9,:),[1,size(FTAvgArray, 4)]);
BispectrumAvg = reshape(BSAvgArray(1,3,9,:,:), [size(BSAvgArray, 4),size(BSAvgArray, 5)]);

eta = sqrt(var(Tau));
eta_true = eta;

estimated_noise_sigma = sqrt(mean([PowerSpectrumAvg(1:N) PowerSpectrumAvg(end-N+1:end)])/(2*N));
disp('True Noise Level:')
disp(true_noise_sigma)
disp('Estimated Noise Level:')
disp(estimated_noise_sigma)

if strcmp(RandomDilationOpts.DilationCalc,'Oracle')
    %eta = eta_true;
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
elseif strcmp(RandomDilationOpts.DilationCalc,'Empirical')
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

%should estimate eta here
[X_w,Y_w] = meshgrid(w); 
if(true_noise_sigma == 0)
    data = 4 .* BispectrumAvg +  CalculateDerivatives(BispectrumAvg, X_w, Y_w, spacing);    
else
    %estimation of noise via power spectrum
    GaussianWidth = GaussianConstant *(estimated_noise_sigma/M).^(1/6);
    unbiased_BS =  UnbiasBispectrum(BispectrumAvg, FourierTransformAvg, estimated_noise_sigma, N);
    %make low pass
    low_pass = MakeSmoothingMatrixBS(w, GaussianWidth);
    low_pass = low_pass/sum(abs(low_pass), 'all');
    %filter with low pass via fft
    unbiased_BS_fft =  fft2(unbiased_BS);
    low_pass_fft =  fft2(low_pass);
    lowpass_BS = ifftshift(ifft2(low_pass_fft.*unbiased_BS_fft));
    
    %d_r(f \ast \phi)
    lowpass_BS_deriv = CalculateDerivatives(lowpass_BS, X_w, Y_w, spacing);
    data = 4 * lowpass_BS + lowpass_BS_deriv; 
    %data = 4 * unbiased_BS + lowpass_BS_deriv; 
end

C_0 = (1-sqrt(3)*eta)/(1+sqrt(3)*eta);
C_1 = 2*sqrt(3)*eta;
C_2 = 1/(1+sqrt(3)*eta);

%True signal
if true_noise_sigma == 0
    real_g0 = real(BispectrumAvg); 
    imag_g0 = imag(BispectrumAvg);
else
    real_g0 = real(lowpass_BS);
    imag_g0 = imag(lowpass_BS);
end
fun_real = @(real_g)compute_loss_and_grad(X_w, Y_w, w, real_g, eta, real(data),RandomDilationOpts,0);
fun_imag = @(imag_g)compute_loss_and_grad(X_w, Y_w, w, imag_g, eta, imag(data),RandomDilationOpts,0);
tol = OptimizationOpts.Uniformtol;

bs_options = optimoptions('fminunc', ...
                       'SpecifyObjectiveGradient',true, ...
                       'MaxFunctionEvaluations', 100000, ...
                       'HessianApproximation', {"lbfgs",100}, ...
                       'MaxIterations',10000, ...
                       'StepTolerance',tol, ...
                       'FunctionTolerance', tol, ...
                       'OptimalityTolerance',tol, ...
                       'Display','iter','StepTolerance', tol);

[est_bispectrum_real,lossval_real] = fminunc(fun_real,real_g0,bs_options);
[est_bispectrum_im,lossval_imag] = fminunc(fun_imag,imag_g0,bs_options);

PlotFourierTransformAvg = FourierTransformAvg .* 2^l;
PlotUnbiasedPS = UnbiasedPS.*(2^(2*l));
PlotTargetPowerSpectrum= 2^(2*l) .* TargetPowerSpectrum;
Plotrecovered_bs = complex(est_bispectrum_real, est_bispectrum_im) .*(2^l).^3;

[f_est_aligned_no_ub_aps, rel_error_inv_no_ub_aps] = invert_bs(PlotFourierTransformAvg, 2^(2*l) .* TargetPowerSpectrum,lowpass_BS.*(2^(3*l)), f, l, 'APS');
[f_est_aligned_aps, rel_error_inv_aps] = invert_bs(PlotFourierTransformAvg, PlotUnbiasedPS, Plotrecovered_bs, f, l,'APS');

gcf = figure
hold on
plot(t, f_est_aligned_no_ub_aps * 2^l)
plot(t, f_est_aligned_aps * 2^l)
plot(t,f)
xlim([-4,4])
no_ub_str = sprintf('NO UB: %.2f', rel_error_inv_no_ub_aps);
ub_str = sprintf('UB: %.2f', rel_error_inv_aps);
legend({no_ub_str, ub_str, 'Ground Truth'}, 'fontsize', 12)
hold off
if(strcmp(RandomDilationOpts.DilationCalc, 'Oracle') && (true_noise_sigma==0.5))
    file_name = sprintf('f%d_inv_plot_oracle_sigma_half.pdf', Signal);
elseif(strcmp(RandomDilationOpts.DilationCalc, 'Oracle') && (true_noise_sigma==1.0))
    file_name = sprintf('f%d_inv_plot_oracle.pdf', Signal);
elseif(strcmp(RandomDilationOpts.DilationCalc, 'Empirical') && (true_noise_sigma==0.5))
    file_name = sprintf('f%d_inv_plot_empirical_sigma_half.pdf', Signal);
elseif(strcmp(RandomDilationOpts.DilationCalc, 'Empirical') && (true_noise_sigma==1.0))
    file_name = sprintf('f%d_inv_plot_empirical.pdf', Signal);
end
exportgraphics(gcf, file_name)


