function Y = CenterBS(Bf)
%This function centers the bispectrum so frequencies range from [-N,N]
%instead of [0,2N]
n = size(Bf,1);
Y = zeros(n,n);
lower = 1:n/2;
upper = n/2+1:n;
Y(upper, upper) = Bf(lower, lower); %top right corner of Y is bottom left corner of Bf
Y(lower, lower) = Bf(upper, upper);
Y(lower, upper) = Bf(upper, lower);
Y(upper, lower) = Bf(lower, upper);
end
