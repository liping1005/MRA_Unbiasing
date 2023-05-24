function [loss, grad] = compute_loss_with_grad_uniform_learn_eta_v2(w,sqrt_g,d,RandomDilationOpts)

g = sqrt_g(1:end-1).^2;
eta = sqrt_g(end);
d_prime = [0 0 diff1_4O_FCD(d,w(2)-w(1)) 0 0];

C_0 = (1-sqrt(3)*eta)/(1+sqrt(3)*eta);
C_1 = 2*sqrt(3)*eta;
C_2 = 1/(1+sqrt(3)*eta);

h = zeros(1,length(g));

if strcmp(RandomDilationOpts.Normalization,'Linf')
    %dilate by C2 and multiply by C1
    h = C_1*(C_2^3)*interp1( w./C_2, d, w, RandomDilationOpts.InterpolationMethod);
    %confused on how to double dilate?
    dilated_h_C_0 = (C_0^3)*interp1( w./C_0, h, w, RandomDilationOpts.InterpolationMethod);
elseif strcmp(RandomDilationOpts.Normalization,'L1')
    %dilate by C2 and multiply by C1
    h = C_1*(C_2)*interp1( w./C_2, d, w, RandomDilationOpts.InterpolationMethod);
    %dilate by C0 now
    dilated_h_C_0 = (C_0)*interp1( w./C_0, h, w, RandomDilationOpts.InterpolationMethod);
end

dilated_d_prime = C_1 * C_2^3 .* interp1( w./C_2, d_prime, w,RandomDilationOpts.InterpolationMethod);
%derivative of L_C0 C1 L_C2 d with respect to the variable z = C0 C2 w
dilated_h_C_0_prime = [0 0 diff1_4O_FCD(dilated_h_C_0, w(2) - w(1)) 0 0];
%corresponds to g = (I-L_C0) C1 L_C2 d
loss = sum((g -h - dilated_h_C_0).^2); 
    
%%% Now compute gradient %%%
    
% First gradient wrt g:
f = 2*(g -h - dilated_h_C_0);

no_interp_idx = find((w<C_0*C_2*min(w)) | (w>C_0*C_2*max(w)));
f(no_interp_idx)=0; %set to zero values which can't be interpolated
% grad_g = 2A^* f where A = I - L_{C_0}
grad_sqrt_g = 2.*sqrt_g(1:end-1).* f;

% Next gradient wrt eta:
if strcmp(RandomDilationOpts.Normalization,'Linf')
    %derivative
    C0_prime = -2*sqrt(3)*(1+sqrt(3)*eta)^(-2);
    C1_prime = 2*sqrt(3);
    C2_prime = -sqrt(3)*(1+sqrt(3)*eta)^(-2);
    %%% individual terms %%%
    %first term: C_1' L_C2 d = C1' h/C1
    %second term: C1 (L_C2 d)' = C1 ((C2)^3 d(C2 w))' = 3 C1 (C2)^3 d(C2 w) + C1 (C2)^3 C2' d'(C2 w) 
    %                                                 = 3h +  LC_2 d' * C2'* w 
    %third term: 3 C0' C0^2 h
    %fourth term: C0^3 dh/dz dz/dn = C0^3 dh/dz *(C2 * C0 w)/dn
    %                              = C0^3 dh/dz C2'C0 w + C0'C2 w
    partial_f_eta = C1_prime*h/C_1;
    partial_f_eta = partial_f_eta + 3*h + (C2_prime.* w) .* dilated_d_prime;
    partial_f_eta = partial_f_eta + 3 * C0_prime * C_0^2  .* dilated_h_C_0/C_0^3;
    partial_f_eta = partial_f_eta + dilated_h_C_0_prime .*((C_2 *C0_prime + C_0 * C2_prime) * w);
elseif strcmp(RandomDilationOpts.Normalization,'L1') %don't handle this case yet
    disp('Not Implemented Error')
    exit   
end
        
grad_eta = sum(f.*partial_f_eta);
     
%Concatenate to get final gradient:
grad = [grad_sqrt_g grad_eta];
    
end