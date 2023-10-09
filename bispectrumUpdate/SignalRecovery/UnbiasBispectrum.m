function UnbiasedBS = UnbiasBispectrum(bispectrum, signal_mean, sigma, N)
    %{
    Calculates estimate for centering of bispectrum based on paper.
    Assumes the bipectrum is a square.
    ARGS
    bispectrum: is the average bispectrum of f_j
    signal_mean: is the mean of the fourier transforms
    sigma: is the noise level
    %}
    [J,~] = size(bispectrum);
    m_x = round(J/2);
    
    %hat_mu(0) part
    mu = signal_mean(m_x+1)/(2*N); %this is equal to mean(NoisyPaddedDilatedSignals(:))
    %noise part
    noise = sigma^2;

    UnbiasingMat = zeros(J);
    UnbiasingMat(1,:) = ones(size(bispectrum(1,:)));
    UnbiasingMat(:,1) = ones(size(bispectrum(:,1)));
    UnbiasingMat = UnbiasingMat + eye(J);
    UnbiasingMat(1,1) = 3; 
    UnbiasingMat = fftshift(UnbiasingMat);
    UnbiasedBS = bispectrum - (2*N)^2 * noise * mu * UnbiasingMat;
end

%avg value where bias is divide by mu * (2N)^2 
%N/2 N/4 N/8 for length?



    
