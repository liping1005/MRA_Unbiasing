function [loss, grad] = compute_loss_with_grad_uniform(w,sqrt_g,C_0,h,RandomDilationOpts)

g = sqrt_g.^2;

%Approximate loss:
if strcmp(RandomDilationOpts.Normalization,'L1')
    dilated_g = C_0*interp1( w./C_0, g, w,RandomDilationOpts.InterpolationMethod);
elseif strcmp(RandomDilationOpts.Normalization,'Linf')
    dilated_g = (C_0^3)*interp1( w./C_0, g, w,RandomDilationOpts.InterpolationMethod);
end
if sum(ismember(w,0))==0
   dilated_g(w<min(w./C_0))=0; %This is only needed for WSC; for PS condition never satisfied
end
loss = sum((g - dilated_g - h).^2);


% Find values where dilation by 1/C_0 is a problem:
no_interp_idx = find((w<C_0*min(w)) | (w>C_0*max(w)));

%Compute gradient:
f = g - dilated_g - h;

if strcmp(RandomDilationOpts.Normalization,'Linf')
    % %Compute gradient. Shrink the functions separately (this does not seem to
    % % help anything)
    %shrink_f1 = ((C_0)^2)*interp1( w.*C_0, g - dilated_g, w,RandomDilationOpts.InterpolationMethod);
    %shrink_f2 = ((C_0)^2)*interp1( w.*C_0, h, w,RandomDilationOpts.InterpolationMethod);
    %shrink_f = shrink_f1-shrink_f2;
    shrink_f = ((C_0)^2)*interp1( w.*C_0, f, w,RandomDilationOpts.InterpolationMethod);

    % %Compute gradient. Avoid stacking interpolations (this does not seem to
    % % help anything)
    %
    % (to use revert to
    % compute_loss_with_grad_uniform(w,sqrt_g,C_0,C_1,C_2,h,d,RandomDilationOpts) )
    %
    % %Find values where dilation by C_2/C_0 is a problem:
    % no_interp_idx2 = find((w<(C_0/C_2)*min(w)) | (w>(C_0/C_2)*max(w)));
    % 
    % %Compute grad term part 1:
    % shrink_g = (C_0^2)*interp1( w.*C_0, g, w,RandomDilationOpts.InterpolationMethod);
    % shrink_g(no_interp_idx) = 0;
    % grad1 = 2*((1+C_0^5)*g - shrink_g - dilated_g);
    % 
    % %Compute grad term part 2:
    % stretch_d = 2*(C_1)*(C_2^3)*interp1( w./C_2, d, w,RandomDilationOpts.InterpolationMethod);
    % shrink_d = 2*(C_0^2)*(C_1)*(C_2^3)*interp1( w.*(C_0/C_2), d, w,RandomDilationOpts.InterpolationMethod);
    % shrink_d(no_interp_idx2) = 0;
    % grad2 = shrink_d - stretch_d;
    % 
    % grad = grad1+grad2;
    
elseif strcmp(RandomDilationOpts.Normalization,'L1')
    shrink_f = interp1( w.*C_0, f, w,RandomDilationOpts.InterpolationMethod);
end

shrink_f(no_interp_idx)=0;
grad = 4*sqrt_g.*(f-shrink_f);


end

