%https://bitbucket.org/annavlittle/inversion-unbiasing/src/master/SupportingFunctions/MakeSmoothingMatrix.m
function SmoothingMatrixBS = MakeSmoothingMatrixBS(w, width)
    SmoothingMatrixBS = zeros(length(w),length(w));
    %width = 10;
    for i=1:length(w)
        for j=1:length(w)
            SmoothingMatrixBS(i,j) = (1/(2*pi*width^2))*exp( (-w(i)^2 - w(j)^2) / (2*width^2) );
        end
    end

