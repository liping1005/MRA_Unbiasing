function [d1] = diff1_4O_FCD(x,h)
%second order centered finite difference method for 4th derivative

%coefficients: 1/12	?2/3	0	2/3	?1/12
n = length(x); %the 4th deriv will have size n-4

%should be able to get derivs for points 3:n-2. I'll need these values:
d1 = ((1/12)*x(1:n-4) + (-2/3)*x(2:n-3) + (0)*x(3:n-2) + (2/3)*x(4:n-1) + (-1/12)*x(5:n))./h;

end
