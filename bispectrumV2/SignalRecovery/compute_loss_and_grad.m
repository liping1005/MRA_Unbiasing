function [loss, grad] = compute_loss_and_grad(X_w, Y_w, w, g, eta, d,RandomDilationOpts, printres)
%https://bitbucket.org/annavlittle/inversion-unbiasing/src/master/SupportingFunctions/compute_loss_with_grad_uniform_learn_eta.m

%{
Inputs
------
w: frequencies sampled
curr_g: current guess for solution to optimization problem
eta: variance of dilation distribution
d: data term, pre transformed?
RandomDilationOpts: dilation operations
%}

%constants from notation section in paper
C_0 = (1-sqrt(3)*eta)/(1+sqrt(3)*eta);
C_1 = 2*sqrt(3)*eta;
C_2 = 1/(1+sqrt(3)*eta);

%Approximate loss:
if strcmp(RandomDilationOpts.Normalization,'L1')
    dilated_g = C_0.*interp2( X_w./C_0, Y_w./C_0, g, X_w, Y_w,RandomDilationOpts.InterpolationMethod);
    dilated_d = C_2.*interp2( w./C_2, w./C_2, d, w, w,RandomDilationOpts.InterpolationMethod);
elseif strcmp(RandomDilationOpts.Normalization,'Linf')
    dilated_g = (C_0^4).*interp2( X_w./C_0, Y_w./C_0, g, X_w, Y_w,RandomDilationOpts.InterpolationMethod);
    dilated_d = (C_2^4).*interp2( X_w./C_2, Y_w./C_2, d, X_w, Y_w,RandomDilationOpts.InterpolationMethod);
end
if printres == 1
    figure
    hold on
    grid on
    axis square
    imagesc('XData',w,'YData',w,'CData',dilated_g)
    title('Real Estimate Dilated Iteration','fontsize',16,'Interpreter','Latex')
    colorbar()
end

%if sum(ismember(w,0))==0 %sanity check basically
%   dilated_g(w<min(w./C_0), :)=0; 
%   dilated_g(:, w<min(w./C_0))=0;
%   dilated_d(w<min(w./C_2),:)=0;
%   dilated_d(:,w<min(w./C_2))=0;
%end

%corresponds to loss in paper
loss = sum((g - dilated_g - C_1 .* dilated_d).^2, 'all');

% Find values where dilation by 1/C_0  or 1/C_2 is a problem:
no_interp_idx = find((w<C_0*min(w)) | (w>C_0*max(w)) | (w>C_2*max(w)) | w<C_2*min(w));

%Compute gradient
%g-dilatated_g corresponds to Ag in paper
%C_1*dilated_d is C_1*L_(C_2)*d from the paper
f = g - dilated_g - C_1.*dilated_d;

%corresponds to C_0^2 h(w/C_0) from adjoint
%recalculate adjoint
if strcmp(RandomDilationOpts.Normalization,'Linf')
    shrink_f = ((C_0)^2).*interp2( C_0.*X_w, C_0.*Y_w, f, X_w, Y_w,RandomDilationOpts.InterpolationMethod);

elseif strcmp(RandomDilationOpts.Normalization,'L1')
    shrink_f = interp2( X_w.*C_0, Y_w.*C_0, f, X_w, Y_w,RandomDilationOpts.InterpolationMethod);
end


shrink_f(no_interp_idx,:)=0;
shrink_f(:,no_interp_idx)=0;
%2(f-shrink_f) is the adjoint/gradient term
grad = 2.*(f-shrink_f);
end