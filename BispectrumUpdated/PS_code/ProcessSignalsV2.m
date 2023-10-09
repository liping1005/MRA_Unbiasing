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
    Undilatedf_hat = ifftshift(fft(fftshift(f))).*(1/2^l);
    UndilatedPowerSpectrum  = abs(Undilatedf_hat).^2;
        
    % Compute Bispectrum for Dilated Signals
    FourierTransformAvg = zeros(1,length(w));
    PowerSpectrumAvg = zeros(1,length(w));
    %disp(size(FourierTransformAvg))
    BispectrumAvg = zeros(length(w), length(w));  
    %disp(size(BispectrumAvg))
    for i=1:M
        ft_temp = ifftshift(fft(fftshift(NoisyPaddedDilatedSignals(i,:)))).*(1/2^l);
        FourierTransformAvg = FourierTransformAvg + ft_temp; 
        PowerSpectrumAvg = PowerSpectrumAvg + abs(ft_temp).^2;
    end  
    FourierTransformAvg = FourierTransformAvg/M;
    PowerSpectrumAvg = PowerSpectrumAvg/M;
elseif strcmp(RandomDilationOpts.SynthesisDomain, 'Frequency')  
    
    w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
    spacing = w(2) - w(1);
    Undilatedf_hat = f1(w);
    %PS
    UndilatedPowerSpectrum  = abs(Undilatedf_hat).^2;
    %BS
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
    PowerSpectrumAvg = mean(abs(FourierTransform).^2);
end