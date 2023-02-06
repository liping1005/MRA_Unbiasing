SmoothingMatrix = zeros(length(w),length(w));
CorrectionMatrix = zeros(length(w),length(w));
DerivSmoothingMatrix = zeros(length(w),length(w));
NormalizationConstant = zeros(1,length(w));
if true_noise_sigma>0
    %width = 3*( (sqrt(2^l)*estimated_noise_sigma)^4/(M*(w(2)-w(1))^2) )^(1/6); %for v5 simulation
    %width = 5*( (sqrt(2^l)*estimated_noise_sigma)^4/(M*(w(2)-w(1))^2) )^(1/6); %for v4, v7 simulation
    %width = ( (sqrt(2^l)*estimated_noise_sigma)^4/(M*(w(2)-w(1))^2) )^(1/6); %for v9-v12 simulation (after renormalizing all signals to have norm 1)
    %width = .2*(M)^(-1/6);
    %width = 10*M^(-1/2);
    %width = M^(-1/3);
    %width = M^(-1/4);
    width = 30*( (sqrt(2^l)*estimated_noise_sigma)^4/(M*(w(2)-w(1))^2) )^(1/3); %setting for simulations v1-v3
else
    %width = 30*( (sqrt(2^l)*0.05)^4/(M*(w(2)-w(1))^2) )^(1/3); %v6 simulation
    %width = 5*( (sqrt(2^l)*0.0625)^4/(M*(w(2)-w(1))^2) )^(1/6); %v8 simulation
end
for i=1:2*N*2^l
    SmoothingMatrix(:,i) = (1/sqrt(2*pi*width^2))*( exp(-(w-w(i)).^2/(2*width^2))+exp(-(w-(2^l)*2*pi-w(i)).^2/(2*width^2))+exp(-(w+(2^l)*2*pi-w(i)).^2/(2*width^2)) )*( w(2)-w(1) );
    NormalizationConstant(i) = sum(SmoothingMatrix(:,i));
    DerivSmoothingMatrix(:,i) = (1/sqrt(2*pi*width^2))*( ((-2*w+2*w(i))/(2*width^2)).*exp(-(w-w(i)).^2/(2*width^2))+((-2*(w-(2^l)*2*pi)+2*w(i))/(2*width^2)).*exp(-(w-(2^l)*2*pi-w(i)).^2/(2*width^2))+((-2*(w+(2^l)*2*pi)+2*w(i))/(2*width^2)).*exp(-(w+(2^l)*2*pi-w(i)).^2/(2*width^2)) )*( w(2)-w(1) );
    CorrectionMatrix(:,i) = SmoothingMatrix(:,i) + (w'-w(i)).*DerivSmoothingMatrix(:,i);
end