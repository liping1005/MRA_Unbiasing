function [ftau] = DilateFunction(f, t, tau)
%Input: 
%    @f: function to dilate
%    tau: M-dimensional column vector of dilation factors (f_tau(x) = f(x(1-\tau))
%    t: N-dimensional row vector of input values at which to evaluate the dilation
%Output:
%   ftau: M by N array where row i is the function dilated by tau(i)
%   evaluated for all values in t
    M = length(tau); %number of dilations
    N = length(t); %number of time points
    ftau = zeros(M,N);
    for i=1:M
        ftau(i,:) = f( (1/(1-tau(i))).*t );
    end
end