% Take Derivatives:
Deltaw = abs(w(2)-w(1));
if strcmp(RandomDilationOpts.SmoothPS, 'yes')
    %FirstDerivPowerSpectrum = [0 0 diff1_4O_FCD(SmoothingMatrix*MeanPowerSpectrum',Deltaw)' 0 0];
    FirstDerivPowerSpectrum = (DerivSmoothingMatrix*(PowerSpectrumAvg)')';
    if strcmp(RandomDilationOpts.SmoothDerivPS, 'yes')
        FirstDerivPowerSpectrum = (SmoothingMatrix*FirstDerivPowerSpectrum')'; %if want to smooth Derivative also
    end
    FirstDerivPowerSpectrum = FirstDerivPowerSpectrum(3:end-2); %to get dimensions to match; maybe recode everything full dimensional
    if PlotFigs==1
        figure
        plot(w(3:end-2),FirstDerivPowerSpectrum)
        title('Derivative of Smoothed PS','Fontsize',14)
    end
elseif strcmp(RandomDilationOpts.SmoothPS, 'no')
    FirstDerivPowerSpectrum = diff1_4O_FCD(PowerSpectrumAvg,Deltaw);
end

% Compute data terms:
if strcmp(RandomDilationOpts.Normalization,'Linf')
    data_term_PS = 3*TargetPowerSpectrum + [0 0 w(3:end-2).*FirstDerivPowerSpectrum 0 0]; 
    if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm, 'yes')
        CorrectionTerm = (CorrectionMatrix*(PowerSpectrumAvg' - 2*N*estimated_noise_sigma^2))';
        data_term_PS_corrected = 3*TargetPowerSpectrum + [0 0 w(3:end-2).*FirstDerivPowerSpectrum 0 0] - CorrectionTerm; 
    end        
elseif strcmp(RandomDilationOpts.Normalization,'L1') 
    data_term_PS = TargetPowerSpectrum + [0 0 w(3:end-2).*FirstDerivPowerSpectrum 0 0]; 
end

%eta_initialization = .35+(.1*rand-0.05); % \eta is initialized Unif(0.3,0.4) 
%eta_init_grid = 0.3;
eta_init_grid = .1:.01:.35;
%eta_initialization = 0.30;    
% Unbias the PS (if differentiable):
if strcmp(RandomDilationOpts.PSUniform,'no')
    UnbiasedPS = TargetPowerSpectrum; 
elseif strcmp(RandomDilationOpts.PSUniform,'yes')            
    % Optimize to find the G = (g, \eta) which produces data_term   
    fun = @(sqrt_g)compute_loss_with_grad_uniform_learn_eta(w,sqrt_g,data_term_PS,RandomDilationOpts);
    tol = OptimizationOpts.Uniformtol;
    %Constrained optimization
    A = zeros(length(TargetPowerSpectrum)+1,length(TargetPowerSpectrum)+1);   
    A(1,length(TargetPowerSpectrum)+1) = 1;
    A(2,length(TargetPowerSpectrum)+1) = -1;
    b = zeros(length(TargetPowerSpectrum)+1,1);
    b(1) = 0.4; %upper bound for eta
    %b(2) = 0;
    b(2) = -.05; %negative of lower bound for eta
    ps_options = optimoptions('fmincon','Algorithm','interior-point', 'SpecifyObjectiveGradient',true,'MaxFunctionEvaluations', 100000,'MaxIterations',10000,'StepTolerance',tol,'FunctionTolerance', tol,'OptimalityTolerance',tol,'Display','iter');
    lossvals_grid_PS = zeros(size(eta_init_grid));
    eta_grid_PS = zeros(size(eta_init_grid));
    UnbiasedPS_grid = zeros(length(eta_init_grid),length(TargetPowerSpectrum));
    for i=1:length(eta_init_grid)
        eta_initialization = eta_init_grid(i);
        sqrt_g0 = [sqrt(abs(TargetPowerSpectrum)) eta_initialization]; 
        %sqrt_g0 = [sqrt(abs(randn(size(TargetPowerSpectrum)))) eta_initialization]; 
        [G,lossval,exitflag,output,lambda,grad,hessian] = fmincon(fun,sqrt_g0,A,b,[],[],[],[],[],ps_options);
        lossvals_grid_PS(i) = lossval;
        eta_grid_PS(i) = G(length(G));
        UnbiasedPS_grid(i,:) = G(1:length(G)-1).^2;
    end
    cand_eta_idx = intersect( find(eta_grid_PS > 0.1) , find(eta_grid_PS < 0.35) ); % only consider interior solutions, throw away eta's too close to the boundary
    if isempty(cand_eta_idx)
        best_eta_idx = find(lossvals_grid_PS==min(lossvals_grid_PS)); % If no interior solutions exist, use the one with smallest lossval
    else
        best_eta_idx = cand_eta_idx( find(lossvals_grid_PS(cand_eta_idx)==min(lossvals_grid_PS(cand_eta_idx))) ); % take one with minimum loss value if there is more than one interior solution
    end
    eta_PS = eta_grid_PS(best_eta_idx);
    UnbiasedPS = UnbiasedPS_grid(best_eta_idx,:);
         
        %%
    if PlotFigs==1
        figure
        plot(w,UnbiasedPS)
        hold on
        plot(w,UndilatedPowerSpectrum)
        plot(w,TargetPowerSpectrum)
        legend({'Unbiased','True','MeanPS'},'fontsize',14)
    end
end