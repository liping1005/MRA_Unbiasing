%https://bitbucket.org/annavlittle/inversion-unbiasing/src/master/SupportingFunctions/ProcessSignals.m

%Randomly sample dilation factors:
Tau = zeros(M,1);
if strcmp(RandomDilationOpts.Distribution,'NoDilation')
    Tau = zeros(M,1);
elseif strcmp(RandomDilationOpts.Distribution,'Uniform')==1
    Tau=2*RandomDilationOpts.MagnitudeMaxTau*rand(M,1)-RandomDilationOpts.MagnitudeMaxTau;
elseif strcmp(RandomDilationOpts.Distribution,'TruncatedGaussian')==1
    pd = makedist('Normal');
    pd.mu = 0;
    pd.sigma = RandomDilationOpts.SDTruncatedGaussian;
    truncGaus = truncate(pd,-RandomDilationOpts.MagnitudeMaxTau,RandomDilationOpts.MagnitudeMaxTau);
    Tau = random(truncGaus,M,1);
end

t1=-(N/2):(1/2^l):(N/2)-1/2^l;
t = -(N):(1/2^l):(N)-1/2^l;

% Dilate Signals, either in space or frequency:
if strcmp(RandomDilationOpts.SynthesisDomain, 'Space')
    % Dilate Signals
    DilatedSignals = zeros(M,length(t1));
    %DilatedSignals = DilateFunction(f1,t1,Tau); %Don't L^1 normalize dilations
    if strcmp(RandomDilationOpts.Distribution,'NoDilation')==1 
        DilatedSignals = ones(M,1)*f1(t1);
    else
        if strcmp(RandomDilationOpts.Normalization,'L1')
            DilatedSignals = (1./(1-Tau)).*DilateFunction(f1,t1,Tau);
        elseif strcmp(RandomDilationOpts.Normalization,'Linf')
            DilatedSignals = DilateFunction(f1,t1,Tau);
        end 
    end

    % Pad with zeros:
    f = [zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)];
    PaddedDilatedSignals = [zeros(M,(2^l)*N/2) DilatedSignals zeros(M,(2^l)*N/2)];
    if strcmp(RandomDilationOpts.Translate, 'True')
        for i=1:M
            rand_trans = randsample(length(t),1);
            PaddedDilatedSignals(i,:) = circshift( PaddedDilatedSignals(i,:), rand_trans );
        end
    end

    % Add Additive Noise
    NoisyPaddedDilatedSignals = PaddedDilatedSignals + true_noise_sigma*sqrt(2^l)*randn( size(PaddedDilatedSignals) );

        
    % Define frequencies, switch to bispectrum here
    w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
    spacing = w(2) - w(1);

    %get fft
    Undilatedf_hat = ifftshift(fft(fftshift(f)));
    UndilatedPowerSpectrum = abs(Undilatedf_hat).^2;
    UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));
        
    % Compute Bispectrum for Dilated Signals
    FourierTransformArr = zeros(M,length(w));
    PowerSpectrumArr = zeros(M,length(w));
    BispectrumAvg = zeros(length(w), length(w));  
    
    for i=1:M
        FourierTransformArr(i,:) = ifftshift(fft(fftshift(NoisyPaddedDilatedSignals(i,:))))*(1/2^l); 
        PowerSpectrumArr(i,:) = abs(FourierTransformArr(i,:)).^2;
        BispectrumAvg = BispectrumAvg + CenterBS(ComputeBispectrum(fftshift(FourierTransformArr(i,:))));
    end 
    FourierTransformAvg = mean(FourierTransformArr);
    PowerSpectrumAvg = mean(PowerSpectrumArr);
    BispectrumAvg = BispectrumAvg/M;
elseif strcmp(RandomDilationOpts.SynthesisDomain, 'Frequency')  
    w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
    spacing = w(2) - w(1);
    Undilatedf_hat = f1(w);
    %PS
    UndilatedPowerSpectrum  = abs(Undilatedf_hat).^2;
    %BS
    UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));
    f = ifftshift(ifft(ifftshift(Undilatedf_hat)));    
    FourierTransform = zeros(M,length(w));
    if strcmp(RandomDilationOpts.Distribution,'NoDilation')==1 
        FourierTransform = ones(M,1)*Undilatedf_hat + true_noise_sigma*randn( size(FourierTransform) );
    else
        if strcmp(RandomDilationOpts.Normalization,'L1')
            FourierTransform = DilateFunction(f1,w,Tau./(Tau-1)) + true_noise_sigma*sqrt(2*N)*randn( size(FourierTransform) );
        elseif strcmp(RandomDilationOpts.Normalization,'Linf')
            FourierTransform = (1-Tau).*DilateFunction(f1,w,Tau./(Tau-1)) + true_noise_sigma*sqrt(2*N)*randn( size(FourierTransform) );
        end 
    end
    %switch to bispectrum here
    BispectrumAvg = zeros(length(w), length(w)); 
    for i=1:M
        BispectrumAvg = BispectrumAvg + CenterBS(ComputeBispectrum(fftshift(FourierTransform(i,:))));
    end 
    BispectrumAvg = BispectrumAvg/M;
    FourierTransformAvg = mean(FourierTransform);
    PowerSpectrumAvg = mean(abs(FourierTransform).^2);
end


eta_true = sqrt(var(Tau));

estimated_noise_sigma = sqrt(mean([PowerSpectrumAvg(1:N) PowerSpectrumAvg(end-N+1:end)])/(2*N));
disp('True Noise Level:')
disp(true_noise_sigma)
disp('Estimated Noise Level:')
disp(estimated_noise_sigma)

if strcmp(RandomDilationOpts.DilationCalc,'Oracle')
    eta = eta_true;
    if strcmp(RandomDilationOpts.SmoothPS, 'yes')
        MakeSmoothingMatrixPS
        TargetPowerSpectrum = (SmoothingMatrix*(PowerSpectrumAvg - 2*N*estimated_noise_sigma^2)')';
    elseif strcmp(RandomDilationOpts.SmoothPS, 'no')
        TargetPowerSpectrum = (PowerSpectrumAvg - 2*N*estimated_noise_sigma^2);
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
        TargetPowerSpectrum =(PowerSpectrumAvg - 2*N*estimated_noise_sigma^2);
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

    %try taking derivative of low pass
    %low_pass_deriv = CalculateDerivatives(low_pass, X_w, Y_w, spacing); 
    %low_pass_deriv = low_pass_deriv/sum(abs(low_pass_deriv), 'all');
    %low_pass_deriv_fft = fft2(low_pass_deriv);
    %lowpass_BS_deriv = ifftshift(ifft2(low_pass_deriv_fft.*unbiased_BS_fft));

    
    %d_r(f \ast \phi), probably too noisy
    lowpass_BS_deriv = CalculateDerivatives(lowpass_BS, X_w, Y_w, spacing);
    %data = 4 * lowpass_BS + lowpass_BS_deriv; 
    data = 4 * lowpass_BS + lowpass_BS_deriv; 
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