userpath('/mnt/ufs18/home-109/yinlipi1/bispectrumV2')
addpath(genpath('../../NumericalExperiments/Oracle'))
addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

plot_settings

f1 = @(x)(4.45458)*(sinc(.2.*(x-32))+sinc(.2.*(-x-32)));
RandomDilationOpts.SynthesisDomain = 'Frequency';

setup

if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm,'yes')
    save('results_f4_oracle.mat','Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','BSerrorVec_corrected','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDBSerrorVec_corrected','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec','SDLogBSerrorVec_corrected')
elseif strcmp(RandomDilationOpts.DilationCalc,'Oracle')
    save('results_f4_oracle.mat','Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec')
elseif strcmp(RandomDilationOpts.DilationCalc,'Empirical')
    save('results_f4_oracle.mat','Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec','Eta_PS_Mean','Eta_PS_Error','Eta_PS_AllSims')  
end