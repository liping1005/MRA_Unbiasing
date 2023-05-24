function deriv = CalculateDerivatives(matrix, X_w, Y_w ,h)
    [len_x,len_y] = size(matrix);
    dy = ((1/12)*matrix(1:len_x-4,:) + (-2/3)*matrix(2:len_x-3,:) + (0)*matrix(3:len_x-2,:) + (2/3)*matrix(4:len_x-1, :) + (-1/12)*matrix(5:len_x,:))./h;
    dx = ((1/12)*matrix(:,1:len_y-4) + (-2/3)*matrix(:,2:len_y-3) + (0)*matrix(:,3:len_y-2) + (2/3)*matrix(:,4:len_y-1) + (-1/12)*matrix(:,5:len_y))./h;
    grad_y = padarray(dy,[2 0], 0);
    grad_x = padarray(dx,[0 2], 0);
    %calculate r d/dr
    %[grad_x, grad_y] = gradient(matrix, h);
    deriv = X_w.* grad_x + Y_w .* grad_y;

