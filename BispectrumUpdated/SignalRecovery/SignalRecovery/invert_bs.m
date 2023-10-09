function [f_est_aligned, rel_error] = invert_bs(fft_guess, unbiasedPS, recoveredBS, f, l, mode) 
    % Now recover f with bispectrum inversion and Boumal's code
    %mean_est = mean(f);
    mean_est = mean(ifftshift(ifft(ifftshift(fft_guess))))./ 2^(l);
    P_est = fftshift(unbiasedPS)'./ 2^(2*l);
    B_est = CenterBS(recoveredBS)./ 2^(3*l);

    fft_shift_avg = fftshift(fft_guess)./ 2^(l);
    if strcmp(mode, 'APS')
        z_est = phases_from_bispectrum_APS_real(B_est, fft_shift_avg(:,1));
    elseif strcmp(mode, 'FM')
        z_est = phases_from_bispectrum_FM_real(B_est);
    elseif strcmp(mode, 'Manifold')
        z_est =  phases_from_bispectrum_complex(B_est);
    elseif strcmp(mode, 'LLL')
        z_est = phases_from_bispectrum_LLL_real(CenterBS(real(recoveredBS)));
    else
        disp('ERROR. Pick Valid Option.')
        exit
    end

    % Recombine (this is computationally cheap)
    f_est = combine_features(mean_est, P_est, z_est);
    f_est = f_est';

    %center 
    f_est_aligned = align_to_reference(f_est, f);
    %rel error
    rel_error = norm(transpose(f)-f_est_aligned*2^l,2)/(norm(transpose(f),2));
end
