BSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
BSerrorVec = zeros( length(Mvalues), 1);
SDBSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
SDBSerrorVec = zeros( length(Mvalues), 1);
SDLogBSerrorVec_NoDilUnbias = zeros( length(Mvalues), 1);
SDLogBSerrorVec = zeros( length(Mvalues), 1);

for s=1:length(Mvalues)
    M = Mvalues(s);
    temp_BSerror_NoDilUnbias = zeros( NumberSimulationsPerValue, 1);
    temp_BSerror = zeros( NumberSimulationsPerValue, 1);
    if strcmp(RandomDilationOpts.DilationCalc,'Empirical')
        temp_eta_BS = zeros( NumberSimulationsPerValue, 1);
        temp_eta_BS_error = zeros( NumberSimulationsPerValue, 1);
    end
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