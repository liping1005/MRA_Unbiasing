function loss = compute_loss_fixed_g(w,sqrt_g,d, eta_grid, RandomDilationOpts)

g = sqrt_g;
h = zeros(1,length(g));
loss = zeros(length(eta_grid), 1);

for i=1:size(eta_grid,2)
    C_0 = (1-sqrt(3)*eta_grid(1,i))/(1+sqrt(3)*eta_grid(1,i));
    C_1 = 2*sqrt(3)*eta_grid(1,i); 
    C_2 = 1/(1+sqrt(3)*eta_grid(1,i));
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
    if (i == 5)
        figure
        plot(w, h)
        hold on
        plot(w, g - dilated_g)
        hold off
        legend({'Dilated Data Term','Dilated Ground Truth'},'fontsize',14) 
    end
    loss(1,i) = sum((g - dilated_g - h).^2); 
end
   
    



    