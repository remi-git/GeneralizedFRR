function [kappa_fit, beta, var_corr, err_k, err_beta] = PSD_fit(trajectory, fminFFT, facq, gamma, kT, begin_fit, end_fit, varargin)
% Computes the Fourier tranform and Lorentzian
% fit of 'trajectory', returns PSD parameters: 
% stiffness and conversion factor
options.MaxIterations = 50;
options.Algorithm = 'levenberg-marquardt';
options.FunctionTolerance = 1e-15;
    
Nfft =  floor(facq/fminFFT);% Nombre of points in the FFT
[spectrum, f] = pwelch(trajectory,  hanning(Nfft), floor(Nfft*0.8), Nfft, facq, 'twosided');
% Model:
fun = @(param,freq)param(1)./(2 * pi^2*(param(2).^2 + freq.^2));
initial_param = [1; 100];
% Fit result:
[P,~,~,covB,MSE] = nlinfit(f(begin_fit:end_fit),...
    spectrum(begin_fit:end_fit), fun, initial_param, options);

P(2) = abs(P(2));
ErrGamma = 0.05;
D_fit = abs(P(1));
kappa_fit = 2 * pi * gamma * abs(P(2));
beta = sqrt( kT / (D_fit * gamma));
var_corr = D_fit/(2*pi*P(2)); % systematic error on variance
err_k = sqrt( (2 * pi * gamma * sqrt(covB(2,2)))^2 + (2 * pi * P(2) * ErrGamma * gamma)^2);
err_Omega0 =  2 * pi * sqrt(covB(2,2));
disp("kappa = " + num2str(1e6*kappa_fit, 4) + " +/- " + num2str(1e6*err_k, 3) + " [pN/micron]");
err_D = sqrt(covB(1,1));
err_beta = sqrt( (sqrt(kT/(gamma)) * D_fit^(-3/2) * err_D * D_fit)^2 + (sqrt(kT/(D_fit)) * gamma^(-3/2) * ErrGamma * gamma)^2);
disp("Fit err D = " + num2str(err_D*100, 4) + " % ");
disp("err gamma = " + num2str(ErrGamma*100, 4) + " % ");
disp("err beta = " + num2str(err_beta/beta*100, 4) + " % ");

for k= 1 : size(varargin,2)
    switch lower(varargin{k})
        case 'verbose'
            figure(1)
            loglog(f, spectrum, 'o');
            hold on;
            loglog(f,P(1)./(2 * pi^2*(P(2).^2 + f.^2)),'b--','LineWidth', 1.3);
            xline(P(2), 'k', 'LineWidth', 2)
            yL = get(gca,'YLim');
            x = [fminFFT begin_fit begin_fit fminFFT];
            y = [yL(1) yL(1) yL(2) yL(2)];
            patch(x,y,'k', 'Facealpha', 0.2)
            x = [end_fit max(f) max(f) end_fit];
            y = [yL(1) yL(1) yL(2) yL(2)];
            patch(x,y,'k', 'Facealpha', 0.05)
            xlabel('Hz')
            ylabel('PSD [V^2/Hz]')
            legend('Data','Lorentzian fit')
            axis tight
            xlim([1, 2e4])
            set(gca,'FontSize',15)
            saveas(gca,'PSD','jpg')
    end
end

