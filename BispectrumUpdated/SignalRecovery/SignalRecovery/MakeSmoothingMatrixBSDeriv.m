%https://bitbucket.org/annavlittle/inversion-unbiasing/src/master/SupportingFunctions/MakeSmoothingMatrix.m
function SmoothingMatrix = MakeSmoothingMatrixBSDeriv(w, width)
    SmoothingMatrix = zeros(length(w),length(w));
    %width = 10;
    for i=1:length(w)
        for j=1:length(w)
            SmoothingMatrix(i,j) = -(1/(4*pi*width^4))*(w(i)^2 + w(j)^2)^(1/2) * exp( (-w(i)^2 - w(j)^2) / (2*width^2) );
        end
    end