function loss = compute_loss(w,sqrt_g,d,RandomDilationOpts,eta)
g = sqrt_g(1:end).^2;
C_0 = (1-sqrt(3)*eta)/(1+sqrt(3)*eta);
C_1 = 2*sqrt(3)*eta;
C_2 = 1/(1+sqrt(3)*eta);


h = zeros(1,length(g));

if strcmp(RandomDilationOpts.Normalization,'Linf')
    
    h = C_1*(C_2^3)*interp1( w./C_2, d, w, RandomDilationOpts.InterpolationMethod);
    dilated_g = (C_0^3)*interp1( w./C_0, g, w,RandomDilationOpts.InterpolationMethod);
    
elseif strcmp(RandomDilationOpts.Normalization,'L1')
    
    h = C_1*(C_2)*interp1( w./C_2, d, w, RandomDilationOpts.InterpolationMethod);
    dilated_g = (C_0)*interp1( w./C_0, g, w,RandomDilationOpts.InterpolationMethod);
    
end
    
% For WSC, set to zero the NAN interpolation values
if sum(ismember(w,0))==0

    no_interp_idx = find(w<min(w./C_2));
    h(no_interp_idx)=0;
    dilated_g(no_interp_idx)=0;

end

loss = sum((g - dilated_g - h).^2); 