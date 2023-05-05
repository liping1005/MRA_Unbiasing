function [loss, grad] = compute_loss_with_grad_uniform_learn_eta(w,sqrt_g,d,RandomDilationOpts)

g = sqrt_g(1:end-1).^2;
eta = sqrt_g(end);
g_prime = [0 0 diff1_4O_FCD(g,w(2)-w(1)) 0 0];
d_prime = [0 0 diff1_4O_FCD(d,w(2)-w(1)) 0 0];

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
   
dilated_d = interp1( w./C_2, d, w,RandomDilationOpts.InterpolationMethod);
dilated_g_prime = interp1( w./C_0, g_prime, w,RandomDilationOpts.InterpolationMethod);
dilated_d_prime = interp1( w./C_2, d_prime, w,RandomDilationOpts.InterpolationMethod);
    
% For WSC, set to zero the NAN interpolation values
if sum(ismember(w,0))==0

    no_interp_idx = find(w<min(w./C_2));
    h(no_interp_idx)=0;
    dilated_g(no_interp_idx)=0;

end

loss = sum((g - dilated_g - h).^2); 
    
% Now compute gradient
    
% First gradient wrt g:
f = g - dilated_g - h;

if strcmp(RandomDilationOpts.Normalization,'Linf')
    shrink_f = ((C_0)^2)*interp1( w.*C_0, f, w,RandomDilationOpts.InterpolationMethod);
elseif strcmp(RandomDilationOpts.Normalization,'L1')
    shrink_f = interp1( w.*C_0, f, w,RandomDilationOpts.InterpolationMethod);
end

no_interp_idx = find((w<C_0*min(w)) | (w>C_0*max(w)));
shrink_f(no_interp_idx)=0; %set to zero values which can't be interpolated
% grad_g = 2A^* f where A = I - L_{C_0}
grad_sqrt_g = 4*sqrt_g(1:end-1).*( f - shrink_f );

% Next gradient wrt eta:
if strcmp(RandomDilationOpts.Normalization,'Linf')
    
    B0 = 2/(1+sqrt(3)*eta)^5;
    B1 = -3*(sqrt(3)-3*eta)*(-1+3*eta^2);
    B2 = sqrt(3)-9*eta*(1-sqrt(3)*eta+eta^2);
    B3 = -sqrt(3)+3*eta+6*sqrt(3)*eta^2;
    B4 = 3*eta;
    partial_f_eta = B0*(B1*dilated_g/(C_0^3)+B2*w.*dilated_g_prime+B3*dilated_d+B4*w.*dilated_d_prime);
    
elseif strcmp(RandomDilationOpts.Normalization,'L1')
    
    B0 = 2/(1+sqrt(3)*eta)^3;
    B1 = sqrt(3)+3*eta;   
    B2 = sqrt(3)-3*eta;
    B3 = -(sqrt(3)+3*eta);
    B4 = 3*eta;   
    partial_f_eta = B0*(B1*dilated_g/(C_0)+B2*w.*dilated_g_prime+B3*dilated_d+B4*w.*dilated_d_prime);
    
end
        
grad_eta = sum(2*f.*partial_f_eta);
     
%Concatenate to get final gradient:
grad = [grad_sqrt_g grad_eta];
    
end