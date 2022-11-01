function deriv = CalculateDerivatives(matrix, X_w, Y_w ,h)
   %calculate r d/dr
   [grad_x, grad_y] = gradient(matrix,h);
   deriv = X_w.* grad_x + Y_w .*grad_y;


