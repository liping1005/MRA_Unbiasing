%https://bitbucket.org/annavlittle/inversion-unbiasing/src/master/SupportingFunctions/MakeSmoothingMatrix.m
function SmoothingMatrix = MakeSmoothingMatrixBS(w)
    SmoothingMatrix = zeros(length(w),length(w));
    %width = 1/10;
    width = 0.5;
    for i=1:length(w)
        for j=1:length(w)
            SmoothingMatrix(i,j) = (1/sqrt(2*pi*width^2))*exp( (-w(i)^2 - w(j)^2)/(2*width^2) );
        end
    end

