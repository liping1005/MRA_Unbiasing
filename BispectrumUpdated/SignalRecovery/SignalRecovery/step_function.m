function [y] = step_function(x,a,b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

y = ones(size(x));
too_big_idx = find(x > b);
too_small_idx = find(x < a);
y(too_big_idx) = zeros(size(too_big_idx));
y(too_small_idx) = zeros(size(too_small_idx));

end