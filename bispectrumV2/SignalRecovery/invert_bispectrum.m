function invert_bispectrum(fft_guess, unbiasedPS, recoveredBS, f, w, t, l)
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [.1 .2 .8 .4]);
    subplot(1,2,1)
    imagesc('XData',w,'YData',w,'CData',real(recoveredBS))
    colorbar
    title('Bispectrum (real part)','FontSize',16)
    set(gca,'YDir','normal')

    % Now recover f with bispectrum inversion and Boumal's code
    %mean_est = mean(f);
    mean_est = real(mean(ifftshift(ifft(ifftshift(fft_guess)))));
    P_est = fftshift(unbiasedPS)';
    B_est = CenterBS(real(recoveredBS));

    fft_shift_avg = fftshift(fft_guess);
    z_est = phases_from_bispectrum_FM_real(B_est, real(fft_shift_avg(:,1)));
    %f_init = combine_features(mean_est, P_est, z_est);

    %f_fft_init = fft(fftshift(f_init))';
    %z_init = exp(1j*angle(f_fft_init));
    %z_est = phases_from_bispectrum_real(B_est, sign(mean_est), z_init');

    % Recombine (this is computationally cheap)
    f_est = combine_features(mean_est, P_est, z_est);
    f_est = f_est';

    disp(size(f_est))

    %center
    f_est_aligned = align_to_reference(f_est, f);
    %rel error
    rel_error = norm(transpose(f)-f_est_aligned*2^l,2)/(norm(f,2) + 1e-6);

    subplot(1,2,2)
    plot(t, f, 'linewidth', 1.5)
    hold on
    plot(t, f_est_aligned.*2^l,'linewidth', 1.5)
    %plot(t_ext,f_est,'*','linewidth',1.5)
    title(sprintf('Rel Error: %.2f', rel_error),'FontSize',12)

