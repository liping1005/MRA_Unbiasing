userpath('/mnt/ufs18/home-109/yinlipi1/bispectrum')
addpath(genpath('../../BS_inversion'))
addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../NumericalExperiments'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

Signal = 14;

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
elseif Signal ==11
    w = gausswin(10, 2.5);
    w = w/sum(w);
    f1 = @(x) 1.2693 .* filter(w,1,sqrt(zigzag((x+10)/5)+zigzag((x-10)/5)));
    RandomDilationOpts.SynthesisDomain = 'Frequency'; 
elseif Signal ==12;
    w = gausswin(10, 2.5);
    w = w/sum(w);
    f1 = @(x) 1.0028 .* filter(w,1,sqrt(zigzag((x+1)/8)+zigzag((x-1)/8)));
    RandomDilationOpts.SynthesisDomain = 'Frequency'; 
elseif Signal == 14
    f1 = @(x)(4.45458)*(sinc(.2.*(x-12))+sinc(.2.*(-x-12))); 
    RandomDilationOpts.SynthesisDomain = 'Frequency';
else 
    disp('Error')
end

plot_settings_learn_eta_sigma_half
setup_learn_eta

if(true_noise_sigma == 0.5)
    myStr = sprintf('results_f%d_l%d_sigma_half_learn_eta.mat', Signal, l);
elseif(true_noise_sigma == 1.0)
    myStr = sprintf('results_f%d_l%d_learn_eta.mat', Signal, l);
end
if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm,'yes')
    save(myStr,'Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','BSerrorVec_corrected','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDBSerrorVec_corrected','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec','SDLogBSerrorVec_corrected', 'SignalerrorVec', 'SDSignalerrorVec',  'SDLogSignalerrorVec')
elseif strcmp(RandomDilationOpts.DilationCalc,'Oracle')
    save(myStr,'Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec', 'SignalerrorVec', 'SDSignalerrorVec',  'SDLogSignalerrorVec')
elseif strcmp(RandomDilationOpts.DilationCalc,'Empirical')
    save(myStr,'Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec', 'SignalerrorVec', 'SDSignalerrorVec',  'SDLogSignalerrorVec', 'Eta_PS_Mean','Eta_PS_Error','Eta_PS_AllSims')
end