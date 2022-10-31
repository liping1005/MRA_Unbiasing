function deriv = CalculateDerivatives(matrix, w)
   %calculate r d/dr
   X = repmat(w, length(w), 1);
   Y = repmat(flip(w)',1,length(w));
   %W_X = repmat(t, length(t), 1);
   %W_Y = repmat(flip(t)',1,length(t));
   %matrix_fft   = fft(fftshift(matrix));
   %dx_fft = 1i * W_X .* matrix_fft;
   %dy_fft = 1i * W_Y .* matrix_fft;
   %dx = ifft(ifftshift(dx_fft));
   %dy = ifft(ifftshift(dy_fft));
   [dx,dy] = gradient(matrix, w(2)-w(1));
   deriv = X .* dx + Y .* dy;


