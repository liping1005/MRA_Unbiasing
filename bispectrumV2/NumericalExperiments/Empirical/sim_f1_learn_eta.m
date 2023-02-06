addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../NumericalExperiments/Empirical'))
addpath(genpath('../../NumericalExperiments/Oracle'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../Experiments'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

plot_settings_learn_eta

% Medium frequency example:
f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);

setup

if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm,'yes')
    save('results_f1_learn_eta.mat','Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','BSerrorVec_corrected','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDBSerrorVec_corrected','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec','SDLogBSerrorVec_corrected')
elseif strcmp(RandomDilationOpts.DilationCalc,'Oracle')
    save('results_f1_learn_eta.mat','Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec')
elseif strcmp(RandomDilationOpts.DilationCalc,'Empirical')
    save('results_f1_learn_eta.mat','Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec','Eta_PS_Mean','Eta_PS_Error','Eta_PS_AllSims')  
end

PlotScriptUniformUnbiasBS
