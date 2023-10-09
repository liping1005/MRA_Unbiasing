addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

Signal = input(['\nWhich example do you want to run?',...
       '\n ',...
       '\n1 Low Low Frequency Gabor',...
       '\n ',...
       '\n2 Low Frequency Gabor',...
       '\n ',...
       '\n3 Medium Frequency Gabor',...
       '\n ',...
       '\n4 Step Function',...
       '\n ',...
       '\n5 Squared Sinc',...
       '\n ',...
       '\n6 Triangle',...
       '\n ',...
       '\n7 Zero Function',...
       '\n ',...
       '\nInput example number without brackets or parentheses: ']);

if Signal==1
    f1 = @(x)(9.759)*exp(-5*x.^2).*cos(4.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';  
elseif Signal==2
    f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==3 % consider making 16 the highest?
    f1 = @(x)(10.6857)*exp(-5*x.^2).*cos(16.*x);
    RandomDilationOpts.SynthesisDomain = 'Space'; 
elseif Signal==4
    f1 = @(x)(4.0 * step_function(x,-1.0,1.0));
    RandomDilationOpts.SynthesisDomain = 'Space';   
elseif Signal==5
    f1 = @(x)19.6116/2*(sinc(4.*x)).^2;
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==6
    f1 = @(x) 3.472*triangle(x,-1.0, 1.0, 2);
     RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==7
    f1 = @(x)zeros(size(x));
    RandomDilationOpts.SynthesisDomain = 'Space';
else 
    disp('Error')
end

plot_settings
setup

myStr = sprintf('results_f%d_oracle.mat', Signal);
if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm,'yes')
    save(myStr,'Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','BSerrorVec_corrected','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDBSerrorVec_corrected','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec','SDLogBSerrorVec_corrected')
elseif strcmp(RandomDilationOpts.DilationCalc,'Oracle')
    save(myStr,'Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec')
elseif strcmp(RandomDilationOpts.DilationCalc,'Empirical')
    save(myStr,'Mvalues','NumberSimulationsPerValue','BSerrorVec_NoDilUnbias','BSerrorVec','SDBSerrorVec_NoDilUnbias','SDBSerrorVec','SDLogBSerrorVec_NoDilUnbias','SDLogBSerrorVec','Eta_PS_Mean','Eta_PS_Error','Eta_PS_AllSims')  
end

PlotScriptUniformUnbiasBS



    