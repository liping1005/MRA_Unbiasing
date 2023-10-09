%this is for BS 
BSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
BSerrorVec = zeros( length(Mvalues), 1);
SDBSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
SDBSerrorVec = zeros( length(Mvalues), 1);
SDLogBSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
SDLogBSerrorVec = zeros( length(Mvalues), 1);

%this is for actual signal recovery APS
SignalerrorVec_NoDilUnbias_APS = zeros( length(Mvalues), 1);
SignalerrorVec_APS = zeros( length(Mvalues), 1);
SDSignalerrorVec_NoDilUnbias_APS = zeros( length(Mvalues), 1);
SDSignalerrorVec_APS = zeros( length(Mvalues), 1);
SDLogSignalerrorVec_NoDilUnbias_APS = zeros( length(Mvalues), 1);
SDLogSignalerrorVec_APS = zeros( length(Mvalues), 1);

%this is for actual signal recovery FM
SignalerrorVec_NoDilUnbias_FM = zeros( length(Mvalues), 1);
SignalerrorVec_FM = zeros( length(Mvalues), 1);
SDSignalerrorVec_NoDilUnbias_FM = zeros( length(Mvalues), 1);
SDSignalerrorVec_FM = zeros( length(Mvalues), 1);
SDLogSignalerrorVec_NoDilUnbias_FM = zeros( length(Mvalues), 1);
SDLogSignalerrorVec_FM = zeros( length(Mvalues), 1);

%this is to keep all the necessary quantities around.
PS_length = size([zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)], 2);
PSCenteredTensor = zeros(1, NumberSimulationsPerValue, length(Mvalues), PS_length);
BSCenteredTensor = zeros(1, NumberSimulationsPerValue, length(Mvalues), PS_length, PS_length);
UnbiasedPSTensor = zeros(1, NumberSimulationsPerValue, length(Mvalues), PS_length);
UnbiasedBSTensor = zeros(1, NumberSimulationsPerValue, length(Mvalues), PS_length, PS_length);

%load data
if(true_noise_sigma == 0.5)
    myStr =  sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/generatedSignals/f%d_signals_all_l%d_sigma_half.mat', Signal, l);
elseif(true_noise_sigma == 1.0)
    myStr =  sprintf('/mnt/ufs18/home-109/yinlipi1/bispectrum/NumericalExperiments/generatedSignals/f%d_signals_all_l%d.mat', Signal, l);   
end 
load(myStr)


if(strcmp(RandomDilationOpts.SynthesisDomain, 'Space'))
    f = [zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)];
    Undilatedf_hat = ifftshift(fft(fftshift(f)));
    UndilatedPowerSpectrum = abs(Undilatedf_hat).^2;
    UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));
elseif(strcmp(RandomDilationOpts.SynthesisDomain, 'Frequency'))
    Undilatedf_hat = f1(w);
    UndilatedPowerSpectrum  = abs(Undilatedf_hat).^2;
    UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));
    f = ifftshift(ifft(ifftshift(Undilatedf_hat)))*2^l;
end

for s=1:length(Mvalues)
    M = Mvalues(s);
    temp_BSerror_NoDilUnbias = zeros( NumberSimulationsPerValue, 1);
    temp_BSerror = zeros( NumberSimulationsPerValue, 1);
    temp_inversionerror_NoDilUnbias_APS = zeros( NumberSimulationsPerValue, 1);
    temp_inversionerror_APS = zeros( NumberSimulationsPerValue, 1);
    temp_inversionerror_NoDilUnbias_FM = zeros( NumberSimulationsPerValue, 1);
    temp_inversionerror_FM = zeros( NumberSimulationsPerValue, 1);    
    if strcmp(RandomDilationOpts.DilationCalc,'Empirical')
        temp_eta_BS = zeros( NumberSimulationsPerValue, 1);
        temp_eta_BS_error = zeros( NumberSimulationsPerValue, 1);
    end


    for q=1:NumberSimulationsPerValue
        %get specific arrays now 
        FourierTransformAvg = reshape(FTAvgArray(1,q,s,:), [1,size(FTAvgArray, 4)]); 
        PowerSpectrumAvg = reshape(PSAvgArray(1,q,s,:),[1,size(FTAvgArray, 4)]);
        BispectrumAvg = reshape(BSAvgArray(1,q,s,:,:), [size(BSAvgArray, 4),size(BSAvgArray, 5)]);
        RunSimulationBS_learn_eta
        temp_BSerror_NoDilUnbias(q) = BSerror_NoDilUnbias_rel;
        temp_BSerror(q) = BSerror_rel;
        temp_inversionerror_NoDilUnbias_APS(q) = rel_error_inv_no_ub_aps;
        temp_inversionerror_APS(q) = rel_error_inv_aps;
        temp_inversionerror_NoDilUnbias_FM(q) = rel_error_inv_no_ub_fm;
        temp_inversionerror_FM(q) = rel_error_inv_fm;
        PSCenteredTensor(1,q,s,:) = 2^(2*l) .* TargetPowerSpectrum;
        BSCenteredTensor(1,q,s,:,:) = unbiased_BS.*(2^(3*l));
        UnbiasedPSTensor(1,q,s,:) = UnbiasedPS;
        UnbiasedBSTensor(1,q,s,:,:) = recovered_bs;
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
    
    SignalerrorVec_NoDilUnbias_APS(s) =  mean(temp_inversionerror_NoDilUnbias_APS);
    SignalerrorVec_APS(s) = mean(temp_inversionerror_APS);
    SDSignalerrorVec_NoDilUnbias_APS(s) = std(temp_inversionerror_NoDilUnbias_APS);
    SDSignalerrorVec_APS(s) = std(temp_inversionerror_APS);
    SDLogSignalerrorVec_NoDilUnbias_APS(s) = std(log2(temp_inversionerror_NoDilUnbias_APS));
    SDLogSignalerrorVec_APS(s) = std(log2(temp_inversionerror_APS));
    
    SignalerrorVec_NoDilUnbias_FM(s) =  mean(temp_inversionerror_NoDilUnbias_FM);
    SignalerrorVec_FM(s) = mean(temp_inversionerror_FM);
    SDSignalerrorVec_NoDilUnbias_FM(s) = std(temp_inversionerror_NoDilUnbias_FM);
    SDSignalerrorVec_FM(s) = std(temp_inversionerror_FM);
    SDLogSignalerrorVec_NoDilUnbias_FM(s) = std(log2(temp_inversionerror_NoDilUnbias_FM));
    SDLogSignalerrorVec_FM(s) = std(log2(temp_inversionerror_FM));

    if strcmp(RandomDilationOpts.DilationCalc,'Empirical')
        Eta_PS_Mean(s) = mean(temp_eta_PS);
        Eta_PS_Error(s) = mean(temp_eta_PS_error);
    end
end