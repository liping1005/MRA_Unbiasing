function UnbiasedBS = UnbiasBispectrum(bispectrum, signal_mean, sigma)
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
    mean = signal_mean(m_x);
    %noise part
    noise = sigma^2;

    UnbiasingMat = zeros(J);
    UnbiasingMat(1,1) = 2; 
    UnbiasingMat(2,:) = ones(size(bispectrum(2,:)));
    UnbiasingMat(:,2) = ones(size(bispectrum(:,2)));
    UnbiasingMat = UnbiasingMat + eye(J);
    UnbiasingMat = fftshift(UnbiasingMat);
    UnbiasedBS = bispectrum - noise * mean * UnbiasingMat;
end



    
