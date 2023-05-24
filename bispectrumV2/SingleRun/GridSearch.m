addpath(genpath('../../BS_Calculation'))
addpath(genpath('../../SignalRecovery'))
addpath(genpath('../../PS_code'))
addpath(genpath('../../Utils'))

N=2^(4);
l_vals = 3:3;
M_vals = 10^4;
BS_width_vals = 5:5:20;
%this needs to be large
PS_width_vals = 10:5:30;
true_noise_sigma_vals = 0.25:0.25:1.0; %Optional: add additive Gaussian noise

Signal = input(['\nWhich example do you want to run?',...
       '\n ',...
       '\n1 Low Frequency Gabor',...
       '\n ',...
       '\n2 Medium Frequency Gabor',...
       '\n ',...
       '\n3 High Frequency Gabor',...
       '\n ',...
       '\n4 Step Function',...
       '\n ',...
       '\n5 C0 Bump',...
       '\n ',...
       '\n6 Sinc',...
       '\n ',...
       '\n7 Squared Sinc',...
       '\n ',...
       '\n8 Gaussian',...
       '\n ',...
       '\n9 Zero Function',...
       '\n ',...
       '\nInput example number without brackets or parentheses: ']);

if Signal==1  
    f1 = @(x)(10.6768)*exp(-5*x.^2).*cos(8.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==2
    f1 = @(x)(10.6857)*exp(-5*x.^2).*cos(16.*x);
    RandomDilationOpts.SynthesisDomain = 'Space'; 
elseif Signal==3 
    f1 = @(x)(10.6857)*exp(-5*x.^2).*cos(32.*x);
    RandomDilationOpts.SynthesisDomain = 'Space';  
elseif Signal==4
    f1 = @(x)(step_function(x,-1,1)/2.0);
    RandomDilationOpts.SynthesisDomain = 'Space';  
elseif Signal==5  %C0 bump function
    f1 = @(x) bumpC0(x,-N/2,N/2,N/2,0);
    RandomDilationOpts.SynthesisDomain = 'Space';    
elseif Signal==6
    f1 = @(x)sinc(pi.*x)/(3.4519);
    RandomDilationOpts.SynthesisDomain = 'Space';   
elseif Signal==7
    f1 = @(x)8.*(sinc(8.*x)).^2;
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==8
    f1 = @(x) x.^2 .* exp(-x.^2/16)/56.7185;
    RandomDilationOpts.SynthesisDomain = 'Space';
elseif Signal==9
    f1 = @(x)zeros(size(x));
    RandomDilationOpts.SynthesisDomain = 'Space';
else 
    disp('Error')
end

RandomDilationOpts.Normalization = 'Linf'; %Options: L1 or Linf normalized dilations
RandomDilationOpts.Translate = 'True'; 
RandomDilationOpts.MagnitudeMaxTau = 0.5; %Needed for both Uniform and TruncatedGaussian (default: 0.2)
RandomDilationOpts.Distribution='Uniform';
RandomDilationOpts.UnbiasingMethod = 'Uniform';
RandomDilationOpts.PSUniform = 'yes'; %for UnbiasingMethod = 'Uniform'; options: yes or no
RandomDilationOpts.SmoothPS = 'yes';
RandomDilationOpts.SmoothDerivPS = 'no';
RandomDilationOpts.SmoothPSCorrectionTerm = 'no'; %Only compute for Oracle! Hasn't been implemented for Empirical
RandomDilationOpts.InterpolationMethod = 'spline';
RandomDilationOpts.DilationCalc='Oracle';
OptimizationOpts.Method = 'Unconstrained';
OptimizationOpts.Initialization = 'MeanPS_NoDilUnbias'; %options: MeanPS_NoDilUnbias, MeanPS_Order2Unbias, MeanPS_Order4Unbias
OptimizationOpts.tol = 1e-10;
OptimizationOpts.Uniformtol = 1e-10; % when UnbiasingMethod = 'Uniform', tolerance for recovering g from h via optimization
PlotFigs = 'no'; %options: 'yes' or 'no

rel_error = zeros(size(l_vals,2), size(M_vals,2), size(BS_width_vals,2), size(PS_width_vals,2), size(true_noise_sigma_vals,2));
for lv=1:size(l_vals,2)
    for Mv=1:size(M_vals,2)
        for Sv=1:size(true_noise_sigma_vals,2)
            l = l_vals(lv);
            M = M_vals(Mv);
            true_noise_sigma = true_noise_sigma_vals(Sv);
            Tau = zeros(M,1);
            if strcmp(RandomDilationOpts.Distribution,'NoDilation')
                Tau = zeros(M,1);
            elseif strcmp(RandomDilationOpts.Distribution,'Uniform')==1
                Tau=2*RandomDilationOpts.MagnitudeMaxTau*rand(M,1)-RandomDilationOpts.MagnitudeMaxTau; 
            elseif strcmp(RandomDilationOpts.Distribution,'TruncatedGaussian')==1
                pd = makedist('Normal');
                pd.mu = 0;
                pd.sigma = RandomDilationOpts.SDTruncatedGaussian;
                truncGaus = truncate(pd,-RandomDilationOpts.MagnitudeMaxTau,RandomDilationOpts.MagnitudeMaxTau);
                Tau = random(truncGaus,M,1);
            end
            t1=-(N/2):(1/2^l):(N/2)-1/2^l;
            t = -(N):(1/2^l):(N)-1/2^l;
            % Dilate Signals, either in space or frequency:
            if strcmp(RandomDilationOpts.SynthesisDomain, 'Space')
                DilatedSignals = zeros(M,length(t1));
            %DilatedSignals = DilateFunction(f1,t1,Tau); %Don't L^1 normalize dilations
                if strcmp(RandomDilationOpts.Distribution,'NoDilation')==1 
                    DilatedSignals = ones(M,1)*f1(t1);
                else
                    if strcmp(RandomDilationOpts.Normalization,'L1')
                        DilatedSignals = (1./(1-Tau)).*DilateFunction(f1,t1,Tau);
                    elseif strcmp(RandomDilationOpts.Normalization,'Linf')
                        DilatedSignals = DilateFunction(f1,t1,Tau);
                    end 
                end

            % Pad with zeros:
                f = [zeros(1,(2^l)*N/2) f1(t1) zeros(1,(2^l)*N/2)];
                PaddedDilatedSignals = [zeros(M,(2^l)*N/2) DilatedSignals zeros(M,(2^l)*N/2)];
                if strcmp(RandomDilationOpts.Translate, 'True')
                    for i=1:M
                        rand_trans = randsample(length(t),1);
                        PaddedDilatedSignals(i,:) = circshift( PaddedDilatedSignals(i,:), rand_trans );
                    end
                end
               % Add Additive Noise
               NoisyPaddedDilatedSignals = PaddedDilatedSignals + true_noise_sigma*sqrt(2^l)*randn( size(PaddedDilatedSignals) );

        
               % Define frequencies, switch to bispectrum here
               w=-pi*(2^l):(pi/N):pi*(2^l)-(pi/N);
               spacing = w(2) - w(1);

               %get fft
               Undilatedf_hat = ifftshift(fft(fftshift(f)));
               UndilatedPowerSpectrum = abs(Undilatedf_hat).^2;
               UndilatedBispectrum = CenterBS(ComputeBispectrum(fftshift(Undilatedf_hat)));
        
               % Compute Bispectrum for Dilated Signals
               FourierTransformArr = zeros(M,length(w));
               PowerSpectrumArr = zeros(M,length(w));
               BispectrumAvg = zeros(length(w), length(w));  
    
               for i=1:M
                   FourierTransformArr(i,:) = ifftshift(fft(fftshift(NoisyPaddedDilatedSignals(i,:))))*(1/2^l); 
                   PowerSpectrumArr(i,:) = abs(FourierTransformArr(i,:)).^2;
                   BispectrumAvg = BispectrumAvg + CenterBS(ComputeBispectrum(fftshift(FourierTransformArr(i,:))));
               end 
               FourierTransformAvg = mean(FourierTransformArr);
               PowerSpectrumAvg = mean(PowerSpectrumArr);
               BispectrumAvg = BispectrumAvg/M;
            end
            for Bv=1:size(BS_width_vals,2)
                 for Pv=1:size(PS_width_vals,2)
                    GaussianConstant = BS_width_vals(Bv);
                    PSWidthConstant = PS_width_vals(Pv);
                    eta_true = sqrt(var(Tau));

                    estimated_noise_sigma = sqrt(mean([PowerSpectrumAvg(1:N) PowerSpectrumAvg(end-N+1:end)])/(2*N));
                    if strcmp(RandomDilationOpts.DilationCalc,'Oracle')
                        eta = eta_true;
                        if strcmp(RandomDilationOpts.SmoothPS, 'yes')
                            MakeSmoothingMatrixPS
                            TargetPowerSpectrum = (SmoothingMatrix*(PowerSpectrumAvg - 2*N*estimated_noise_sigma^2)')';
                        elseif strcmp(RandomDilationOpts.SmoothPS, 'no')
                            TargetPowerSpectrum = (PowerSpectrumAvg - 2*N*estimated_noise_sigma^2);
                        end
                        if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
                            Unbias_Uniform
                        elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
                            Unbias_GeneralDerivative
                        end
                    elseif strcmp(RandomDilationOpts.DilationCalc,'Empirical')
                        if strcmp(RandomDilationOpts.SmoothPS, 'yes')
                            MakeSmoothingMatrixPS
                            TargetPowerSpectrum = (SmoothingMatrix*(PowerSpectrumAvg - 2*N*estimated_noise_sigma^2)')';
                        elseif strcmp(RandomDilationOpts.SmoothPS, 'no')
                            TargetPowerSpectrum =(PowerSpectrumAvg - 2*N*estimated_noise_sigma^2);
                        end
                        if strcmp(RandomDilationOpts.UnbiasingMethod, 'Uniform') == 1
                            Unbias_Uniform
                        elseif strcmp(RandomDilationOpts.UnbiasingMethod, 'GeneralDerivative') == 1
                            Unbias_GeneralDerivative
                        end  
                    end

                    %should estimate eta here
                    [X_w,Y_w] = meshgrid(w); 
                    if(true_noise_sigma == 0)
                        data = 4 .* BispectrumAvg +  CalculateDerivatives(BispectrumAvg, X_w, Y_w, spacing);    
                    else
                        %estimation of noise via power spectrum
                        GaussianWidth = GaussianConstant *(estimated_noise_sigma/M).^(1/6);
                        unbiased_BS =  UnbiasBispectrum(BispectrumAvg, FourierTransformAvg, estimated_noise_sigma, N);
                        %make low pass
                        low_pass = MakeSmoothingMatrixBS(w, GaussianWidth);
                        low_pass = low_pass/sum(abs(low_pass), 'all');
                        %filter with low pass via fft
                        unbiased_BS_fft =  fft2(unbiased_BS);
                        low_pass_fft =  fft2(low_pass);
                        lowpass_BS = ifftshift(ifft2(low_pass_fft.*unbiased_BS_fft));
                        lowpass_BS_deriv = CalculateDerivatives(lowpass_BS, X_w, Y_w, spacing);
                        %data = 4 * lowpass_BS + lowpass_BS_deriv; 
                        data = 4 * unbiased_BS + lowpass_BS_deriv; 
                    end

                    C_0 = (1-sqrt(3)*eta)/(1+sqrt(3)*eta);
                    C_1 = 2*sqrt(3)*eta;
                    C_2 = 1/(1+sqrt(3)*eta);
                    if true_noise_sigma == 0
                        real_g0 = real(BispectrumAvg); 
                        imag_g0 = imag(BispectrumAvg);
                    else
                        real_g0 = real(lowpass_BS);
                        imag_g0 = imag(lowpass_BS);
                    end
                    fun_real = @(real_g)compute_loss_and_grad(X_w, Y_w, w, real_g, eta, real(data),RandomDilationOpts,0);
                    fun_imag = @(imag_g)compute_loss_and_grad(X_w, Y_w, w, imag_g, eta, imag(data),RandomDilationOpts,0);
                    tol = OptimizationOpts.Uniformtol;

                    bs_options = optimoptions('fminunc', ...
                       'SpecifyObjectiveGradient',true, ...
                       'MaxFunctionEvaluations', 100000, ...
                       'HessianApproximation', {"lbfgs",100}, ...
                       'MaxIterations',10000, ...
                       'StepTolerance',tol, ...
                       'FunctionTolerance', tol, ...
                       'OptimalityTolerance',tol, ...
                       'Display','iter','StepTolerance', tol);

                    [est_bispectrum_real,lossval_real] = fminunc(fun_real,real_g0,bs_options);
                    [est_bispectrum_im,lossval_imag] = fminunc(fun_imag,imag_g0,bs_options);
 
                    PlotFourierTransformAvg = FourierTransformAvg .* 2^l;
                    PlotUnbiasedPS = UnbiasedPS.*(2^(2*l));
                    PlotTargetPowerSpectrum= 2^(2*l) .* TargetPowerSpectrum;
                    Plotrecovered_bs = complex(est_bispectrum_real, est_bispectrum_im) .*(2^l).^3;


                    figure('visible','off');
                    plot(w,PlotUnbiasedPS)
                    hold on
                    plot(w,UndilatedPowerSpectrum)
                    plot(w,PlotTargetPowerSpectrum)
                    legend({'Unbiased','True','MeanPS'},'fontsize',14)
                    hold off
                    
                    fName_PS = sprintf('PS_f%d_l%d_M%d_BSW%d_PSW%d_sigma_%f.fig', Signal, l, int16(log10(M)), GaussianConstant, PSWidthConstant, int16(true_noise_sigma*100));
                    savefig(fName_PS)

                    PlotResults

                    fName_BS = sprintf('BS_f%d_l%d_M%d_BSW%d_PSW%d_sigma_%f.fig', Signal, l, int16(log10(M)), GaussianConstant, PSWidthConstant, int16(true_noise_sigma*100));
                    savefig(fName_BS)
                    
                    inversion
                    fName_inv = sprintf('inv_f%d_l%d_M%d_BSW%d_PSW%d_sigma_%f.fig', Signal, l, int16(log10(M)), GaussianConstant, PSWidthConstant, int16(true_noise_sigma*100));
                    savefig(fName_inv)
                    
                    rel_error(lv, Mv, Bv, Pv, Sv) = norm(f-f_est_aligned)/norm(f);
                end
            end
        end
    end
end
search_str = sprintf('grid_search_f%d.mat', Signal);
save(search_str, 'rel_error');