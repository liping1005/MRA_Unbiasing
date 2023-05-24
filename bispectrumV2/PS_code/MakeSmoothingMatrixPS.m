SmoothingMatrix = zeros(length(w),length(w));
CorrectionMatrix = zeros(length(w),length(w));
DerivSmoothingMatrix = zeros(length(w),length(w));
NormalizationConstant = zeros(1,length(w));
if true_noise_sigma>0
    width =  PSWidthConstant * (estimated_noise_sigma^4/M)^(1/6);
end
for i=1:2*N*2^l
    SmoothingMatrix(:,i) = (1/sqrt(2*pi*width^2))*( exp(-(w-w(i)).^2/(2*width^2))+exp(-(w-(2^l)*2*pi-w(i)).^2/(2*width^2))+exp(-(w+(2^l)*2*pi-w(i)).^2/(2*width^2)) )*( w(2)-w(1) );
    NormalizationConstant(i) = sum(SmoothingMatrix(:,i));
    DerivSmoothingMatrix(:,i) = (1/sqrt(2*pi*width^2))*( ((-2*w+2*w(i))/(2*width^2)).*exp(-(w-w(i)).^2/(2*width^2))+((-2*(w-(2^l)*2*pi)+2*w(i))/(2*width^2)).*exp(-(w-(2^l)*2*pi-w(i)).^2/(2*width^2))+((-2*(w+(2^l)*2*pi)+2*w(i))/(2*width^2)).*exp(-(w+(2^l)*2*pi-w(i)).^2/(2*width^2)) )*( w(2)-w(1) );
    CorrectionMatrix(:,i) = SmoothingMatrix(:,i) + (w'-w(i)).*DerivSmoothingMatrix(:,i);
end