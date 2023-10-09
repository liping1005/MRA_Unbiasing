function [y] = zigzag(x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
y = zeros(size(x));

for i=1:length(x)
    
    if (-3<=x(i)) && (x(i)<-2)
        y(i) = 3 + x(i);
    elseif (-2<=x(i)) && (x(i)<-1)
        y(i) = -1 -x(i);
    elseif (-1<=x(i)) && (x(i)<0)
        y(i) = 1 + x(i);
    elseif (0<=x(i)) && (x(i)<1)
        y(i) = 1 - x(i);
    elseif (1<=x(i)) && (x(i)<2)
        y(i) = -1 + x(i);
    elseif (2<=x(i)) && (x(i)<3)
        y(i) = 3 - x(i);
    else
        y(i) = 0; 
    end
    
end