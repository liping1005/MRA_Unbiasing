function out = bump(x)
    out = zeros(1,length(x));
    for i=1:1:length(x)
        if(abs(x(i)) < 1)
            out(i) = exp(-1/(1-x(i)^2));
        end
    end
end