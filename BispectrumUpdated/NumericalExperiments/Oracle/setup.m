%this is for BS 
BSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
BSerrorVec = zeros( length(Mvalues), 1);
SDBSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
SDBSerrorVec = zeros( length(Mvalues), 1);
SDLogBSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
SDLogBSerrorVec = zeros( length(Mvalues), 1);

%this is for actual signal recovery
SignalerrorVec = zeros( length(Mvalues), 1);
SDSignalerrorVec = zeros( length(Mvalues), 1);
SDLogSignalerrorVec = zeros( length(Mvalues), 1);

%load data
myStr =  sprintf('f%d_signals_all.mat', Signal);
load(myStr)

f = [zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)];
Undilatedf_hat = ifftshift(fft(fftshift(f)));
UndilatedPowerSpectrum = abs(Undilatedf_hat).^2;
UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));

for s=1:length(Mvalues)
    M = Mvalues(s);
    temp_BSerror_NoDilUnbias = zeros( NumberSimulationsPerValue, 1);
    temp_BSerror = zeros( NumberSimulationsPerValue, 1);
    if strcmp(RandomDilationOpts.DilationCalc,'Empirical')
        temp_eta_BS = zeros( NumberSimulationsPerValue, 1);
        temp_eta_BS_error = zeros( NumberSimulationsPerValue, 1);
    end
    %get specific arrays now 
    FourierTransformAvg = reshape(FTAvgArray(1,s,:), [1,size(FTAvgArray, 3)]); 
    PowerSpectrumAvg = reshape(PSAvgArray(1,s,:),[1,size(FTAvgArray, 3)]);
    BispectrumAvg = reshape(BSAvgArray(1,s,:,:), [size(FTAvgArray, 3),size(FTAvgArray, 3)]);

    for q=1:NumberSimulationsPerValue
        RunSimulationBS
        temp_BSerror_NoDilUnbias(q) = BSerror_NoDilUnbias_rel;
        temp_BSerror(q) = BSerror_rel;
        if strcmp(RandomDilationOpts.DilationCalc,'Empirical')
            temp_eta_PS(q) = eta_PS;
            eta = sqrt(var(Tau));
            temp_eta_PS_error(q) = norm(eta_PS-eta);
            Eta_PS_AllSims(s,q) = eta_PS;
        end
    end
    BSerrorVec_NoDilUnbias(s) = mean(temp_BSerror_NoDilUnbias);
    BSerrorVec(s) = mean(temp_BSerror);
    SDBSerrorVec_NoDilUnbias(s) = std(temp_BSerror_NoDilUnbias);
    SDBSerrorVec(s) = std(temp_BSerror);
    SDLogBSerrorVec_NoDilUnbias(s) = std(log2(temp_BSerror_NoDilUnbias));
    SDLogBSerrorVec(s) = std(log2(temp_BSerror));
    if strcmp(RandomDilationOpts.DilationCalc,'Empirical')
        Eta_PS_Mean(s) = mean(temp_eta_PS);
        Eta_PS_Error(s) = mean(temp_eta_PS_error);
    end
end