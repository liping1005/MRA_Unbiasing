ProcessSignalsV2
estimated_noise_sigma = sqrt(mean([PowerSpectrumAvg(1:N) PowerSpectrumAvg(end-N+1:end)])/(2*N));
if strcmp(RandomDilationOpts.SmoothPS, 'yes')
    MakeSmoothingMatrixPS
    TargetPowerSpectrum = (SmoothingMatrix*(PowerSpectrumAvg - 2*N*estimated_noise_sigma^2)')';
elseif strcmp(RandomDilationOpts.SmoothPS, 'no')
    TargetPowerSpectrum = MeanPowerSpectrum - 2*N*estimated_noise_sigma^2;
end


if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
    Unbias_Uniform
elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
    Unbias_GeneralDerivative
end

% Try to recover signal from Avg PS (only works for real,even signal with positive FT)
recoveredf_UnbiasedPS = ifftshift(ifft(fftshift(sqrt(UnbiasedPS))))*2^l;
% Output main quantities of interest:

if norm(UndilatedPowerSpectrum) > 0
    
    PSerror_NoDilUnbias_rel = norm(UndilatedPowerSpectrum - TargetPowerSpectrum)/norm(UndilatedPowerSpectrum);
    if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
        PSerror_rel = norm(UndilatedPowerSpectrum - UnbiasedPS)/norm(UndilatedPowerSpectrum);
        if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm, 'yes')
            PSerror_rel_corrected = norm(UndilatedPowerSpectrum - UnbiasedPS_corrected)/norm(UndilatedPowerSpectrum);
        end
    elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
        PSerror_rel = norm(UndilatedPowerSpectrum(4:end-3) - UnbiasedPS)/norm(UndilatedPowerSpectrum(4:end-3));
    end
    
elseif norm(UndilatedPowerSpectrum) == 0 %for the zero signal, just compute absolute error
    
    PSerror_NoDilUnbias_rel = norm(UndilatedPowerSpectrum - TargetPowerSpectrum);
    if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
        PSerror_rel = norm(UndilatedPowerSpectrum - UnbiasedPS);
        if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm, 'yes')
            PSerror_rel_corrected = norm(UndilatedPowerSpectrum - UnbiasedPS_corrected);
        end
    elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
        PSerror_rel = norm(UndilatedPowerSpectrum(4:end-3) - UnbiasedPS);
    end    
end


%%
%if strcmp(PlotFigs,'yes')==1
    
%    PlotResultsPS
    
%end