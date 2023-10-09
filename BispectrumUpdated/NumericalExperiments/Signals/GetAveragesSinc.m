userpath('/mnt/ufs18/home-109/yinlipi1/bispectrum')
addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

%% Select example: 


Signal = 15;

if Signal==1
    f1 = @(x)(9.759)*exp(-5*x.^2) .* cos(4.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';  
elseif Signal==2
    f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==3 % consider making 16 the highest?
    f1 = @(x)1/(sqrt(0.0035)) .* exp(-32.*x.^2).*cos(16.*x);
    RandomDilationOpts.SynthesisDomain = 'Space'; 
elseif Signal==4 %step function is probably too hard
    f1 = @(x)(4.0 * step_function(x,-1.0,1.0));
    RandomDilationOpts.SynthesisDomain = 'Space';    
elseif Signal==5
    f1 =@(x)1/sqrt(0.0078).*(sinc(4.*x));
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==6
    f1 = @(x) 3.472*triangle(x,-1.0, 1.0, 2);
     RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal == 15
    f1 = @(x) (8.1930* cos(6*x).*step_function(x, -.5, .5));
    RandomDilationOpts.SynthesisDomain = 'Space';
else 
    disp('Error')
end

%% Set Parameters  

N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=4;
M = 2^20;
numMVals = size(4:2:20, 2);
true_noise_sigma_vals = 0.5:1:0.5;
NumSimulations = 5;
RandomDilationOpts.Normalization = 'Linf'; %Options: L1 or Linf normalized dilations
RandomDilationOpts.Translate = 'True'; 
RandomDilationOpts.MagnitudeMaxTau = 0.5; %Needed for both Uniform and TruncatedGaussian (default: 0.2)
RandomDilationOpts.Distribution='Uniform';



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


%%% Dilate The Signals. This stays the same throughout. %%%

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

% Define frequencies
w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
spacing = w(2) - w(1);

% Pad with zeros:
f = [zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)];
PaddedDilatedSignals = [zeros(M,(2^l)*N/2) DilatedSignals zeros(M,(2^l)*N/2)];
if strcmp(RandomDilationOpts.Translate, 'True')
    for i=1:M
        rand_trans = randsample(length(t),1);
        PaddedDilatedSignals(i,:) = circshift( PaddedDilatedSignals(i,:), rand_trans );
    end
end

%get fft
Undilatedf_hat = ifftshift(fft(fftshift(f)));
UndilatedPowerSpectrum = abs(Undilatedf_hat).^2;
UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));



FTAvgArray = zeros(size(true_noise_sigma_vals,2), NumSimulations, numMVals,length(w));
PSAvgArray =  zeros(size(true_noise_sigma_vals,2), NumSimulations, numMVals,length(w));
BSAvgArray =  zeros(size(true_noise_sigma_vals,2), NumSimulations, numMVals,length(w),length(w));

for i=1:NumSimulations
    for j=1:size(true_noise_sigma_vals,2)
        NoisyPaddedDilatedSignals = PaddedDilatedSignals + true_noise_sigma_vals(:,j)*sqrt(2^l)*randn( size(PaddedDilatedSignals) );
        FourierTransformArr = zeros(M,length(w));
        PowerSpectrumArr = zeros(M,length(w));
        BispectrumAvg = zeros(length(w), length(w));  
        %count for which M is saved
        count = 1;
        for k=1:M
            FourierTransformArr(k,:) = ifftshift(fft(fftshift(NoisyPaddedDilatedSignals(k,:))))*(1/2^l); 
            PowerSpectrumArr(k,:) = abs(FourierTransformArr(k,:)).^2;
            BispectrumAvg = BispectrumAvg + CenterBS(ComputeBispectrum(fftshift(FourierTransformArr(k,:))));
            if(ceil(log2(k)) == floor(log2(k)) && (log2(k) >= 4) && mod(floor(log2(k)),2) == 0)
                FourierTransformAvgk = mean(FourierTransformArr(1:k,:));
                PowerSpectrumAvgk = mean(PowerSpectrumArr(1:k,:));
                BispectrumAvgk = BispectrumAvg/k;
             
                %save averages in a file
                FTAvgArray(j,i,count, :) = FourierTransformAvgk;
                PSAvgArray(j,i,count, :) = PowerSpectrumAvgk;
                BSAvgArray(j,i,count, :, :) =  BispectrumAvgk;
                count = count + 1;
            end
        end
    end 
end
strName =  sprintf('f%d_signals_all_l%d_sigma_half.mat', Signal, l);
save(strName, 'FTAvgArray', 'PSAvgArray', 'BSAvgArray','f','Tau')