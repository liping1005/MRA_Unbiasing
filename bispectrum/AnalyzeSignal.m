%Based on https://bitbucket.org/annavlittle/inversion-unbiasing/src/master/SupportingFunctions/AnalyzeSignal.m

function [w, Bispectrum, f_hat] = AnalyzeSignal(inputform,f,N,l)
%Input:
% f: signal defined on [-N, N) with spacing 1/2^l (f has length 2^l*2N; f=f1+noise where f1 lives on [-N/2,N/2))


% Define frequencies:
w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);

% Compute the bispectrum:

if strcmp(inputform,'Signal')
    % Compute fft of signal (note must shift [-N,N) to [0,2N) first):
    % Note: Use Delta t = 1/2^l instead of Delta t = 1 in computation of FT
    f_hat = ifftshift(fft(fftshift(f)))*(1/2^l);
    Bispectrum = CenterBS(ComputeBispectrum(f_hat));
elseif strcmp(inputform,'Bispectrum')
    Bispectrum = f;
end
end