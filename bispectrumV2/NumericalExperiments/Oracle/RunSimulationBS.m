ProcessSignals

%use PS to check if signal is 0
if norm(UndilatedPowerSpectrum) > 0
    BSerror_NoDilUnbias_rel = norm(UndilatedBispectrum - BispectrumAvg)/norm(UndilatedBispectrum);
    if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
        BSerror_rel = norm(UndilatedBispectrum - recovered_bs)/norm(UndilatedBispectrum);
    elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
        BSerror_rel = norm(UndilatedBispectrum - recovered_bs)/norm(UndilatedBispectrum);
    end
    
elseif norm(UndilatedPowerSpectrum) == 0 %for the zero signal, just compute absolute error
    
    BSerror_NoDilUnbias_rel = norm(UndilatedBispectrum - BispectrumAvg);
    if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
        BSerror_rel = norm(UndilatedBispectrum - Unbiased_BS);
    elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
        BSerror_rel = norm(UndilatedBispectrum - Unbiased_BS);
    end    
end