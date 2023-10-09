function [Bf] = ComputeBispectrum(f_hat)
%Input: FT of a signal
%Output: Bispectrum of the signal (according to Singer definition)
Bf1 = (f_hat.')*conj(f_hat);
Bf2 = Circulant(f_hat);
Bf = Bf1.*Bf2;
end

