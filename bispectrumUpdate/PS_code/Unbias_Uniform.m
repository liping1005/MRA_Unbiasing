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
    
%     %For debugging (exact data_term_PS for medium frequency Gabor)
%     g_eta = zeros(size(UndilatedPowerSpectrum));
%     %g_eta(N*2^l+1) = mean((1-Tau).^2)*UndilatedPowerSpectrum(N*2^l+1);
%     for i=1:length(w)
%         if w(i)>0
%             tau = linspace((1-sqrt(3)*sqrt(var(Tau)))*w(i),(1+sqrt(3)*sqrt(var(Tau)))*w(i),1024);
%         elseif w(i)<0
%             tau = linspace((1+sqrt(3)*sqrt(var(Tau)))*w(i),(1-sqrt(3)*sqrt(var(Tau)))*w(i),1024);
%         end
%         g_tau = 2*pi*(exp(-(tau-16).^2./20)+exp(-(tau+16).^2./20)).^2./40;
%         g_eta(i) = sum((tau.^2).*g_tau./((2*sqrt(3*var(Tau)))*abs(w(i))^3))*(tau(2)-tau(1));
%     end
%     g_eta(N*2^l+1) = mean((1-Tau).^2)*UndilatedPowerSpectrum(N*2^l+1);
%     data_term_PS_true = 3*g_eta + [0 0 w(3:end-2).*diff1_4O_FCD(g_eta,Deltaw) 0 0];  
elseif strcmp(RandomDilationOpts.Normalization,'L1') 
    data_term_PS = TargetPowerSpectrum + [0 0 w(3:end-2).*FirstDerivPowerSpectrum 0 0]; 
end


%%
eta_true = sqrt(moment(Tau,2));

if strcmp(RandomDilationOpts.MomentCalc,'Oracle')

    %eta_true = RandomDilationOpts.MagnitudeMaxTau/sqrt(3);
    C_0 = (1-sqrt(3)*eta_true)/(1+sqrt(3)*eta_true);
    C_1 = 2*sqrt(3)*eta_true;
    C_2 = 1/(1+sqrt(3)*eta_true);

    % Unbias the PS (if differentiable):
    if strcmp(RandomDilationOpts.PSUniform,'no')
        UnbiasedPS = TargetPowerSpectrum; %TargetPowerSpectrum = MeanPowerSpectrum - 2*N*noise_sigma^2
    elseif strcmp(RandomDilationOpts.PSUniform,'yes')
        h_PS = zeros(1,length(TargetPowerSpectrum));
        if strcmp(RandomDilationOpts.Normalization,'L1')
            h_PS = C_1*C_2*interp1( w./C_2, data_term_PS, w, RandomDilationOpts.InterpolationMethod);
        elseif strcmp(RandomDilationOpts.Normalization,'Linf')
            h_PS = C_1*(C_2^3)*interp1( w./C_2, data_term_PS, w, RandomDilationOpts.InterpolationMethod);
            if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm, 'yes')
                h_PS_corrected = C_1*(C_2^3)*interp1( w./C_2, data_term_PS_corrected, w, RandomDilationOpts.InterpolationMethod);
            end
        end

        % Optimize to find the g which produces h
        sqrt_g0 = sqrt(abs(TargetPowerSpectrum));
        fun = @(sqrt_g)compute_loss_with_grad_uniform(w,sqrt_g,C_0,h_PS,RandomDilationOpts);
        tol = OptimizationOpts.Uniformtol;
        options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations', 50000,'MaxIterations',5000,'StepTolerance', tol,'FunctionTolerance', tol,'OptimalityTolerance',tol,'Display','iter','CheckGradients',false);
        [sqrt_UnbiasedPS,lossval,exitflag,output,grad]=fminunc(fun,sqrt_g0,options);
        UnbiasedPS = sqrt_UnbiasedPS.^2;
        
        if strcmp(RandomDilationOpts.SmoothPSCorrectionTerm, 'yes')
            sqrt_g0 = sqrt(abs(TargetPowerSpectrum));
            fun = @(sqrt_g)compute_loss_with_grad_uniform_v3(w,sqrt_g,C_0,h_PS_corrected,RandomDilationOpts);
            tol = OptimizationOpts.Uniformtol;
            options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations', 50000,'MaxIterations',5000,'StepTolerance', tol,'FunctionTolerance', tol,'OptimalityTolerance',tol,'Display','iter','CheckGradients',false);
            [sqrt_UnbiasedPS_corrected,lossval,exitflag,output,grad]=fminunc(fun,sqrt_g0,options);
            UnbiasedPS_corrected = sqrt_UnbiasedPS_corrected.^2;
        end
         
        if PlotFigs==1
            figure
            plot(w,UnbiasedPS,'LineWidth',2)
            hold on
            plot(w,UndilatedPowerSpectrum,'LineWidth',2)
            plot(w,TargetPowerSpectrum,'LineWidth',2)
            %plot(w,inverted_WFT,'--','LineWidth',2)
            legend({'Unbiased','True','MeanPS'},'fontsize',14) 
        end


    end   
elseif strcmp(RandomDilationOpts.MomentCalc,'Empirical')
    
    %eta_initialization = .35+(.1*rand-0.05); % \eta is initialized Unif(0.3,0.4) 
    %eta_init_grid = 0.3;
    eta_init_grid = .1:.05:.35;
    %eta_initialization = 0.30;
    
    % Unbias the PS (if differentiable):
    if strcmp(RandomDilationOpts.PSUniform,'no')
        UnbiasedPS = TargetPowerSpectrum; %TargetPowerSpectrum = MeanPowerSpectrum - 2*N*noise_sigma^2
    elseif strcmp(RandomDilationOpts.PSUniform,'yes')
            
       % Optimize to find the G = (g, \eta) which produces data_term   
        fun = @(sqrt_g)compute_loss_with_grad_uniform_learn_eta_v3(w,sqrt_g,sqrt(UndilatedPowerSpectrum),RandomDilationOpts);
        tol = OptimizationOpts.Uniformtol;
        %Constrained optimization
        A = zeros(length(TargetPowerSpectrum)+1,length(TargetPowerSpectrum)+1);   
        A(1,length(TargetPowerSpectrum)+1) = 1;
        A(2,length(TargetPowerSpectrum)+1) = -1;
        b = zeros(length(TargetPowerSpectrum)+1,1);
        b(1) = 0.4; %upper bound for eta
        %b(2) = 0;
        b(2) = -.05; %negative of lower bound for eta
        options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxFunctionEvaluations', 100000,'MaxIterations',10000,'StepTolerance',tol,'FunctionTolerance', tol,'OptimalityTolerance',tol,'Display','off');
        lossvals_grid_PS = zeros(size(eta_init_grid));
        eta_grid_PS = zeros(size(eta_init_grid));
        UnbiasedPS_grid = zeros(length(eta_init_grid),length(TargetPowerSpectrum));
        for i=1:length(eta_init_grid)
            eta_initialization = eta_init_grid(i);
            sqrt_g0 = [sqrt(abs(TargetPowerSpectrum)) eta_initialization]; 
            %sqrt_g0 = [sqrt(abs(randn(size(TargetPowerSpectrum)))) eta_initialization]; 
            [G,lossval,exitflag,output,lambda,grad,hessian] = fmincon(fun,sqrt_g0,A,b,[],[],[],[],[],options);
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
    end    
end
    

%Compute errors:
PSerror = norm(UndilatedPowerSpectrum - UnbiasedPS)/sqrt(2^l);
if and(strcmp(RandomDilationOpts.SmoothPSCorrectionTerm, 'yes'), strcmp(RandomDilationOpts.MomentCalc, 'Oracle'))
   PSerror_corrected = norm(UndilatedPowerSpectrum - UnbiasedPS_corrected)/sqrt(2^l); 
end

%% Plot

if strcmp(PlotFigs,'yes')==1
 
    figure
    
    if strcmp(GlobalOpts.ComputeWavelets,'yes')
        NumFigs = 6;
    else
        NumFigs = 4;
    end

    subplot(NumFigs/2,2,1);
    plot(w,UndilatedPowerSpectrum);
    title('Power Spectrum of Signal','fontsize',14)

    subplot(NumFigs/2,2,2);
    plot(w,PowerSpectrum(1,:));
    title('Power Spectrum of One Noisy Signal','fontsize',14)

    subplot(NumFigs/2,2,3);
    hold on
    for i=1:min(M,20)
        plot(w, PowerSpectrum(i,:))
    end
    xlim([min(w) max(w)])
    xlabel('Frequency \omega','fontsize',14);
    title(['Power Spectrum for M= ' num2str(M) ' Noisy Signals'],'fontsize',14)
    plot(w, UndilatedPowerSpectrum,'Color','k','LineWidth',3)
    
    subplot(NumFigs/2,2,4);
    plot(w, TargetPowerSpectrum,'LineWidth',2)
    hold on
    plot(w,UnbiasedPS,'LineWidth',2)
    plot(w, UndilatedPowerSpectrum,'LineWidth',2)
    title({'Target PS, Mean PS, and Unbiased PS'},'fontsize',14)
    legend({'Mean PS', 'Unbiased PS','Target PS'},'FontSize',14)
    legend('Location','northeast')

    %title({'Undilated Power Spectrum and Average of Dilated Power Spectra', '(With Additive Noise Unbiasing)'},'fontsize',14)

    pos = get(gcf,'position');
    set(gcf,'position',[pos(1:2)/4 pos(3:4)*2])
end