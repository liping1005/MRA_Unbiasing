userpath('/mnt/ufs18/home-109/yinlipi1/bispectrum')
addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

%% Select example: 


Signal = 9;

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
    f1 =@(x)2.0520.*(sinc(x./8));
    RandomDilationOpts.SynthesisDomain = 'Frequency';  
elseif Signal==5
    f1 =@(x)1/sqrt(0.0078).*(sinc(4.*x));
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==6
    f1 = @(x) 3.472*triangle(x,-1.0, 1.0, 2);
     RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==7
    f1 = @(x)zeros(size(x));
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal ==8
    f1=@(x)13.8564 * triangle(x,-0.25, 0, 0.25);
    RandomDilationOpts.SynthesisDomain = 'Space'; 
elseif Signal ==9
    f1=@(x)  6.6109.*(bumpC2(x,-5,-5,2,0) + bumpC2(x,5,5,2,0));
    RandomDilationOpts.SynthesisDomain = 'Frequency'; 
elseif Signal ==11
    w = gausswin(10, 2.5);
    w = w/sum(w);
    f1 = @(x) 1.2693 .* filter(w,1,sqrt(zigzag((x+10)/5)+zigzag((x-10)/5)));
    RandomDilationOpts.SynthesisDomain = 'Frequency'; 
else 
    disp('Error')
end

%% Set Parameters  

N=2^(4); %Choose N at least 8, or we don't get J>0; choose N a power of 2, or weird things happen
l=4;
M = 2^20;
numMVals = size(4:2:20, 2);
true_noise_sigma=1.0;
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

w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
spacing = w(2) - w(1);
Undilatedf_hat = f1(w);
%PS
UndilatedPowerSpectrum  = abs(Undilatedf_hat).^2;
%BS
UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));
f = ifftshift(ifft(ifftshift(Undilatedf_hat)))*2^l; 



FTAvgArray = zeros(1, NumSimulations, numMVals,length(w));
PSAvgArray =  zeros(1, NumSimulations, numMVals,length(w));
BSAvgArray =  zeros(1, NumSimulations, numMVals,length(w),length(w));

for i=1:NumSimulations
        FourierTransform = zeros(M,length(w));
        %create noise in space so it is real valued
        noise = randn(size(FourierTransform));
        noise_fft = fftshift(fft(noise')',2) *(1/2^l);
    
        if strcmp(RandomDilationOpts.Distribution,'NoDilation')==1 
            FourierTransform = ones(M,1)*Undilatedf_hat + true_noise_sigma*noise_fft;
        else
            if strcmp(RandomDilationOpts.Normalization,'L1')
                FourierTransform = DilateFunction(f1,w,Tau./(Tau-1)) + true_noise_sigma*sqrt(2^l)*noise_fft;
            elseif strcmp(RandomDilationOpts.Normalization,'Linf')
                FourierTransform = (1-Tau).*DilateFunction(f1,w,Tau./(Tau-1)) + true_noise_sigma*sqrt(2^l)*noise_fft;
            end 
        end
        BispectrumAvg = zeros(length(w), length(w));  
        %count for which M is saved
        count = 1;
        for k=1:M
            BispectrumAvg = BispectrumAvg + CenterBS(ComputeBispectrum(fftshift(FourierTransform(k,:))));
            if(ceil(log2(k)) == floor(log2(k)) && (log2(k) >= 4) && mod(floor(log2(k)),2) == 0)
                FourierTransformAvgk = mean(FourierTransform(1:k,:));
                PowerSpectrumAvgk = mean(abs(FourierTransform(1:k,:)).^2);
                BispectrumAvgk = BispectrumAvg/k;
                %save averages in a file
                FTAvgArray(1,i,count, :) = FourierTransformAvgk;
                PSAvgArray(1,i,count, :) = PowerSpectrumAvgk;
                BSAvgArray(1,i,count, :, :) =  BispectrumAvgk;
                count = count + 1;
            end
        end 
end

if(true_noise_sigma == 0.5)
    strName =  sprintf('f%d_signals_all_l%d_sigma_half.mat', Signal, l);
elseif(true_noise_sigma == 1.0)
    strName =  sprintf('f%d_signals_all_l%d.mat', Signal, l);
end

save(strName, 'FTAvgArray', 'PSAvgArray', 'BSAvgArray','f','Tau')