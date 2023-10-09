function y = triangle(x, left, right, height)
    N = size(x,2);
    y = zeros(size(x));
    for i=1:N
        if(x(i)<left || x(i)>right)
            y(i) = 0;
        elseif(left < x(i) && x(i) < (left+right)/2)
            slope = height/abs(left);
            y(i) = slope * x(i) + height;
        elseif(right >= x(i)  && x(i) >= (left+right)/2)
            slope = -height/right;
            y(i) = slope * x(i) + height;
        end
    end
end