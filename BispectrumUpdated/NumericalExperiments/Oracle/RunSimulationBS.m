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
%recovered_bs = complex(est_bispectrum_real, est_bispectrum_im);
UnbiasedPS = UnbiasedPS.*(2^l).^2;
recovered_bs = complex(est_bispectrum_real, est_bispectrum_im) .*(2^l).^3;

%figure
%plot(w,UnbiasedPS)
%hold on
%plot(w,UndilatedPowerSpectrum)
%plot(w,2^(2*l) *TargetPowerSpectrum)
%legend({'Unbiased','True','MeanPS'},'fontsize',14)
%hold off

%use PS to check if signal is 0
if norm(UndilatedPowerSpectrum) > 0
    BSerror_NoDilUnbias_rel = norm(UndilatedBispectrum - BispectrumAvg .*2^(3*l))/norm(UndilatedBispectrum);
    if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
        BSerror_rel = norm(UndilatedBispectrum - recovered_bs)/norm(UndilatedBispectrum);
    elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
        BSerror_rel = norm(UndilatedBispectrum - recovered_bs)/norm(UndilatedBispectrum);
    end
    
elseif norm(UndilatedPowerSpectrum) == 0 %for the zero signal, just compute absolute error
    
    BSerror_NoDilUnbias_rel = norm(UndilatedBispectrum - BispectrumAvg .*2^(3*l));
    if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
        BSerror_rel = norm(UndilatedBispectrum - recovered_bs);
    elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
        BSerror_rel = norm(UndilatedBispectrum - recovered_bs);
    end    
end

%invert bispectrum now
%inversion
