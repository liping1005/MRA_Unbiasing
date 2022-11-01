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
    Undilatedf_hat = fft(fftshift(f)).*(1/2^l);
    %compute bispectrum and center
    UndilatedBispectrum = CenterBS(ComputeBispectrum(Undilatedf_hat));
        
    % Compute Bispectrum for Dilated Signals
    FourierTransformAvg = zeros(1,length(w));
    %disp(size(FourierTransformAvg))
    BispectrumAvg = zeros(length(w), length(w));  
    %disp(size(BispectrumAvg))
    for i=1:M
        ft_temp = fft(fftshift(NoisyPaddedDilatedSignals(i,:))).*(1/2^l);
        FourierTransformAvg = FourierTransformAvg + ft_temp; 
        BispectrumAvg = BispectrumAvg + CenterBS(ComputeBispectrum(ft_temp));
    end  
    FourierTransformAvg = FourierTransformAvg/M;
    BispectrumAvg = BispectrumAvg/M;
elseif strcmp(RandomDilationOpts.SynthesisDomain, 'Frequency')  
    
    w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
    spacing = w(2) - w(1);
    Undilatedf_hat = f1(w);
    %switch to bispectrum here
    UndilatedBipectrum = CenterBS(ComputeBispectrum(Undilatedf_hat));
    f = ifftshift(ifft(ifftshift(Undilatedf_hat)))*2^l;    
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
    Bispectrum = CenterBS(ComputeBispectrum(FourierTransform)); 
end

%make guassian filter
low_pass = MakeSmoothingMatrixBS(w);


% Empirically estimate the additive noise variance
% Use somewhere near the corner where the signal would be nearly zero
% basically
estimated_sigma = sqrt(mean(abs([BispectrumAvg(2,1:N) BispectrumAvg(2,end-N+1:end)])));
%estimated_sigma = 0.01;

%calculate \tilde{g}_eta + \tilde{g}_\sigma
unbiased_BS = UnbiasBispectrum(BispectrumAvg, fftshift(FourierTransformAvg), estimated_sigma);

%filter with low pass via fft
unbiased_BS_hat = fft(fftshift(unbiased_BS)).*(1/2^l);
lowpass_hat = fft(fftshift(low_pass)).*(1/2^l);

lowpass_BS = ifft(ifftshift(unbiased_BS_hat .* lowpass_hat)).*2^l;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
%calcuate part of the estimator
%data = 4 * lowpass_BS + r_grid .* CalculateDerivatives(lowpass_BS, X, Y, spacing);

X_w = repmat(w, length(w), 1);
Y_w= repmat(w',1,length(w));
data = 4 .* BispectrumAvg +  CalculateDerivatives(BispectrumAvg, X_w, Y_w,  w(2)-w(1));
%solve optimization problem via whatever optimizer for estimator, here is
%problem
eta = sqrt(var(Tau));
%eta = 0.2;
C_0 = (1-sqrt(3)*eta)/(1+sqrt(3)*eta);
disp(C_0^4)
C_1 = 2*sqrt(3)*eta;
C_2 = 1/(1+sqrt(3)*eta);

%True signal
TrueBS = CenterBS(ComputeBispectrum(Undilatedf_hat));

dilated_g = (C_0^4).*interp2(X_w./C_0, Y_w./C_0, real(TrueBS), X_w, Y_w,RandomDilationOpts.InterpolationMethod);
dilated_d = (C_2^4).*interp2( X_w./C_2, Y_w./C_2, real(data), X_w, Y_w,RandomDilationOpts.InterpolationMethod);


figure
hold on
grid on
axis square
imagesc('XData',w,'YData',w,'CData',real(TrueBS)-dilated_g)
title('I-LC_0','fontsize',16,'Interpreter','Latex')
colorbar()

figure
hold on
grid on
axis square
imagesc('XData',w,'YData',w,'CData',C_1.*dilated_d)
title('C_1 L_C_2','fontsize',16,'Interpreter','Latex')
colorbar()

%what should the initialization be when eta is known?
real_g0 = real(BispectrumAvg); 
imag_g0 = imag(BispectrumAvg);
fun_real = @(real_g)compute_loss_and_grad(X_w, Y_w, w, real_g, eta, real(data),RandomDilationOpts,0);
fun_imag = @(imag_g)compute_loss_and_grad(X_w, Y_w, w, imag_g, eta, imag(data),RandomDilationOpts,0);
tol = OptimizationOpts.Uniformtol;

options = optimoptions('fminunc', ...
                       'SpecifyObjectiveGradient',true, ...
                       'MaxFunctionEvaluations', 100000, ...
                       'HessianApproximation', {"lbfgs",50}, ...
                       'MaxIterations',10000, ...
                       'StepTolerance',tol, ...
                       'FunctionTolerance', tol, ...
                       'OptimalityTolerance',tol, ...
                       'Display','iter','StepTolerance', tol);

[est_bispectrum_real,lossval_real] = fminunc(fun_real,real(TrueBS),options);
[est_bispectrum_im,lossval_imag] = fminunc(fun_imag,imag_g0,options);

%fun_real = @(real_g)obj(w, real_g, eta, real(data),RandomDilationOpts);
%est_bispectrum_real= fminsearch(fun_real,real_g0);

[I,J] = size(TrueBS);
disp(min(real(BispectrumAvg),[],'all'))
disp(min(real(est_bispectrum_real),[],'all'))
disp(min(real(TrueBS),[],'all'))

disp(max(real(BispectrumAvg),[],'all'))
disp(max(real(est_bispectrum_real),[],'all'))
disp(max(real(TrueBS),[],'all'))

disp(sum((real(TrueBS)-real(BispectrumAvg)).^2,'all')/(I*J))
disp(sum((est_bispectrum_real-real(TrueBS)).^2,'all')/(I*J))
disp(sum((real(TrueBS)-real(data)).^2,'all')/(I*J))
%error
%disp( sum((est_bispectrum_real - real(TrueBS)).^2 , 'all')/ sum(real(TrueBS).^2, 'all') )
%CenteredTrueBS = UnbiasBispectrum(TrueBS, f, estimated_sigma);
%use bispectrum inversion algorithm to see if it works