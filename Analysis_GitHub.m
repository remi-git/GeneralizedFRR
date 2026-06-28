%% July 2026
clear all;
clc
disp(' ============================================ ')
disp('This codes loads experimental data ; computes response and correlation functions')
disp('in order to test the equilibrium and generalized FRR - following R. Goerlich et al. Phys Rev E 2026')
disp(' ============================================ ')
disp('')
set(groot, 'defaultFigureCloseRequestFcn', 'close(gcf)');
set(0, 'DefaultFigureRenderer', 'painters');
load('Colors.mat');
%% load steady-state trajectories
disp('Loading Steady-State Data')
Files = [5, 6, 7, 8, 9, 10];
Nfiles = numel(Files);
RawSteadyTraj = zeros(Nfiles,  720000);
Noise = zeros(Nfiles,  720000);
Wb = waitbar(0, 'loading');
for i=1:Nfiles
    waitbar(i/Nfiles)
    % load trajectories calibrated with the python code SteadyStateAnalysis.py
    RawSteadyTraj(i,:) = load(Date + "_RawSteadyTraj_"+Files(i)+".out");
    Noise(i,:) = load(Date + "_Noise_"+Files(i)+".out");
end
delete(Wb)
disp('Done')
%% load calibration data
disp(' ============================================ ')
disp('Stiffness Calibration')
K_calib_1 = load('Values3/Kappa');
V_calib_1 = load('Values3/Voltage');
K_calib_2 = load('Values4/Kappa');
V_calib_2 = load('Values4/Voltage');
MeanKappaCalib = 0.5 * (K_calib_1 + K_calib_2);
Param = polyfit(V_calib_1, MeanKappaCalib, 1);
VCalibArray = linspace(0, 1, 100);

figure()
p1 = plot(V_calib_1, 1e6*K_calib_1, 'o');
p1.MarkerSize = 12;
p1.MarkerFaceColor = Blue;
p1.MarkerEdgeColor = 'k';
hold on
p1 = plot(V_calib_2, 1e6*K_calib_2, 'o');
p1.MarkerSize = 10;
p1.MarkerFaceColor = LightBlue;
p1.MarkerEdgeColor = 'k';
plot(VCalibArray, 1e6*(Param(1) * VCalibArray + Param(2)), 'k--', 'LineWidth', 2)
axis tight
xlim([0.25, 0.75])
ylim([5, 25])
ylabel('$\kappa \rm{ ~ [pN/\mu m]}$', 'interpreter', 'latex')
xlabel('$V_{\rm AOM}  ~\rm{[V]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 20, 'fontname', 'times')
set(gcf, 'units', 'centimeters', 'position', [1 1 12 8]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "Calib.png",'Resolution', 200)
%% Equilibrium traj, measure Kappa, effective temperature
disp('Measure steady-state stiffness and effective temperature')
kB = 1.3806e-23;
T = 293;
eta = 2.414e-5 * 10^(247.5/(T - 140));
R = 1500e-9;
gamma = 6 * pi * eta * R;
D = kB * T / gamma;
facq = 2^15;
dt = 1/facq;
fminFFT = 1;
NCaseErrK = 10;
StartPSD = linspace(2, 5, NCaseErrK);
EndPSD = linspace(1000, 4000, NCaseErrK);
All_Ki_PSD = zeros(NCaseErrK, 1);
All_ErrKi_PSD = zeros(NCaseErrK, 1);
All_BetaFitPsd_Eq = zeros(NCaseErrK, 1);
for i = 1:NCaseErrK % loop over various fit boundaries
    [All_Ki_PSD(i), All_BetaFitPsd_Eq(i), VarCorrEq, All_ErrKi_PSD(i), ErrorBeta] =...
        PSD_fit(RawSteadyTraj(1,:), fminFFT, facq, gamma, kB*T, floor(StartPSD(i)), floor(EndPSD(i)));
end
Ki_PSD = mean(All_Ki_PSD);
ErrKi_PSD = sqrt(var(All_Ki_PSD, 0, 1) ./ NCaseErrK + mean(All_ErrKi_PSD).^2);
BetaFitPsd_Eq = mean(All_BetaFitPsd_Eq);
Si_PSD = BetaFitPsd_Eq^2 * var(RawSteadyTraj(1,:));
Si_PSD_th = kB * T / Ki_PSD;

StartPSD = linspace(5, 10, NCaseErrK);
EndPSD = linspace(1000, 4000, NCaseErrK);
All_Kf_PSD = zeros(NCaseErrK, 1);
All_ErrKf_PSD = zeros(NCaseErrK, 1);
for i = 1:NCaseErrK % loop over various fit boundaries
    [All_Kf_PSD(i), BetaFitPsd, VarCorr, All_ErrKf_PSD(i), ErrorBeta] =...
        PSD_fit(RawSteadyTraj(2,:), fminFFT, facq, gamma, kB*T, floor(StartPSD(i)), floor(EndPSD(i)));
end
Kf_PSD = mean(All_Kf_PSD);
ErrKf_PSD = sqrt( var(All_Kf_PSD, 0, 1) ./ NCaseErrK + mean(All_ErrKf_PSD)^2);
Sf_PSD = BetaFitPsd_Eq^2 * var(RawSteadyTraj(2,:));
Sf_PSD_th = kB * T / Kf_PSD;
SteadyTraj = BetaFitPsd_Eq * RawSteadyTraj; % calibrate trajecotories
disp(' ============================================ ')
%% Load full ensemble equilibrium
disp('Loading Dynamical Data - White Noise')
FilesDyn = [15, 16]; % both experiments
EnsembleEq = load("Values" + FilesDyn(2) + "/" + Date + "_EnsDown_" + FilesDyn(2) +".out");
NumberTrajectories = 18181; %numel(EnsembleEq(1,:));
disp('Done')
FilteredEnsemble = zeros(size(EnsembleEq));
c = 1;
Skip = 100;
disp('Filtering spurious signal at 1757Hz')
Wb = waitbar(0, 'Filtering');
for i = 1:NumberTrajectories
    if i == c
        waitbar(i/NumberTrajectories)
        c = c + Skip;
    end
     % filter out peak at 1757
    FilteredEnsemble(:,i) = NotchFiltering(EnsembleEq(:,i), facq, 1757, 50);
end
delete(Wb)
VarEq = BetaFitPsd_Eq^2 * var(EnsembleEq, 0, 2);
Si_eq = mean(VarEq(80:200));
ErrSi_eq = sqrt(var(VarEq(80:200)));
Sf_eq = mean(VarEq(350:500));
ErrSf_eq = sqrt(var(VarEq(350:500)));
Corr_i = Si_PSD_th / Si_eq;
Corr_f = Sf_PSD_th / Sf_eq;
Corr_var_Eq = (Corr_i + Corr_f) / 2;
Ki_Equipart = kB*T / Si_eq % equipartition-based evaluation of stiffness
ErrKi_equipart = kB*T * ErrSi_eq/Si_eq^2;
Kf_Equipart = kB*T / Sf_eq
ErrKf_Equipart = kB*T * ErrSf_eq/Sf_eq^2;
disp(' ============================================ ')
%% load ensemble colored noise
disp('Loading Dynamical Data - Colored Noise')
EnsembleCn = load("Values" + FilesDyn(1) + "/" + Date + "_EnsDown_" + FilesDyn(1) +".out");
EnsembleNoiseCn = load("Values" + FilesDyn(1) + "/" + Date + "_EnsNoiseDown_" + FilesDyn(1) +".out");
% following 2 lines remove artefacts in the signal: a slight change in mean
% position is observed and substracted
EnsembleNoiseCn(1:200, :) = EnsembleNoiseCn(1:200, :) - mean(mean(EnsembleNoiseCn(1:200, :))); % remove change in mean
EnsembleNoiseCn(201:end, :) = EnsembleNoiseCn(201:end, :) - mean(mean(EnsembleNoiseCn(201:end, :)));
EnsembleNoiseCn(200:202,:) = EnsembleNoiseCn(197:199, :); % remove points lying exaclty at the quench
disp('Done')
FilteredEnsembleCn = zeros(size(EnsembleCn));
c = 1;
Skip = 100;
disp('Filtering spurious signal at 1757Hz')
Wb = waitbar(0, 'Filtering');
for i = 1:numel(EnsembleCn(1,:)) % filter out peak at 1757
    if i == c
        waitbar(i/NumberTrajectories)
        c = c + Skip;
    end
    FilteredEnsembleCn(:,i) = NotchFiltering(EnsembleCn(:,i), facq, 1757, 10);
end
delete(Wb)

TauC = 3e-3; % experimentally imposed noise correlation
OmegaC = 1./TauC;
Omega_i = Ki_Equipart / gamma;
Tau_i = 1/Omega_i;
Omega_f = Kf_Equipart / gamma;
Tau_f = 1/Omega_f;
VarCn = BetaFitPsd_Eq^2 * var(EnsembleCn, 0, 2);
Si_Cn = mean(VarCn(1:200));
Sf_Cn = mean(VarCn(400:end));
SigmaEta_i = sqrt((1 * Si_Cn - D / Omega_i) .* Omega_i * gamma^2 .* (Omega_i + OmegaC));
D2 = (Corr_i * Si_Cn - D * Tau_i) .* (gamma^2 * (TauC + Tau_i) / (TauC^2 * Tau_i^2));
SigmaEta_f = sqrt((1 * Sf_Cn - D / Omega_f) .* Omega_f * gamma^2 .* (Omega_f + OmegaC));
SigmaEta = mean([SigmaEta_i, SigmaEta_f], 2); % Mean noise amplitude (s
disp(' ============================================ ')
%% The end of the time-dependent case is at steady-state : use to evalute Cxx
disp('Evaluate stationnary correlation function')
Xeq_s = BetaFitPsd_Eq * FilteredEnsemble(400:end,:);
X_s = BetaFitPsd_Eq * FilteredEnsembleCn(400:end,:);
Eta_s = SigmaEta * EnsembleNoiseCn(400:end,:) ./ mean(sqrt(var(EnsembleNoiseCn, 0, 2))); % normalized
%% Compute analytical variance
disp('Compute analytical variance')
Ntime = numel(VarEq);
t = 200; % time of step
TimeFull = linspace(-t*dt, (Ntime-t)*dt, Ntime)';
TimeDecay = linspace(0, (Ntime - t)*dt, Ntime - t)';
VarUpAnalytical = zeros(2, Ntime);
% for colored noise, following R. Goerlich et al. Phys Rev E (2026)
InitialNeqVariance = D / Omega_i + SigmaEta^2 / (Omega_i * gamma^2 * (Omega_i + OmegaC));
AnalyticalDecay = InitialNeqVariance .* exp(-2 * Omega_f * TimeDecay)...
        + D / Omega_f .* (1 - exp(-2 * Omega_f * TimeDecay))...
        + SigmaEta^2/gamma^2 * ( 1 / (Omega_f * (Omega_f + OmegaC))...
            + exp(-2 * Omega_f * TimeDecay) / (Omega_f * (Omega_f - OmegaC))...
            - 2 * exp(-(Omega_f + OmegaC) * TimeDecay) / (Omega_f^2 - OmegaC^2) )...
        + 2 * SigmaEta^2/gamma^2 .* (exp(- (Omega_f + OmegaC) * TimeDecay) - exp(-2 * Omega_f * TimeDecay))...
            / ((Omega_i + OmegaC) * (Omega_f - OmegaC));

VarUpAnalytical(1,:) = cat(1, InitialNeqVariance * ones(t, 1), AnalyticalDecay);

% at equilibrium, standard solution of the Fokker-Planck equation
VarUpAnalytical(2,:) = cat(1, D/Omega_i * ones(t, 1),...
    (D/Omega_i - D/Omega_f) .* exp(-2 * Omega_f * TimeDecay) + D/Omega_f);
% Statistical error by Chi2 test
[error_top, error_bottom ] = var_error(VarEq, NumberTrajectories, numel(VarEq), 99.7);
Chi2Error_Eq = (error_top + error_bottom) / 2;
[error_top, error_bottom ] = var_error(VarCn, NumberTrajectories, numel(VarCn), 99.7);
Chi2Error_Cn = (error_top + error_bottom) / 2;
% Plot variance with analytics
Time = linspace(-t * dt, (Ntime-t) * dt, Ntime);
Sp1 = subplot(2, 1, 1);
Sp1.Position = [0.2 0.57 0.75 0.32]; % [left bottom width height]
s = shadedErrorBar(1e3*Time(1:end-1), 1 * VarEq(1:end-1), Chi2Error_Eq,...
    'transparent', true, 'patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle','-', 'Color', Blue)
s.mainLine.LineWidth = 2;
s.mainLine.Color = Blue;
s.patch.FaceColor = Blue;
%errorbar(Time(1:end-1), Corr_var_Eq * VarEq(1:end-1), Chi2Error_Eq, 'o')
hold on
plot(1e3*Time(1:end-1), VarUpAnalytical(2, 1:end-1), 'k--', LineWidth = 3)
axis tight
xlim([-2, 11])
box on
ylabel('$\langle x^2(t) \rangle \rm{ ~ [m^2]}$', 'interpreter', 'latex')
xticklabels('')
set(gca, 'fontsize', 20, 'fontname', 'times')
Sp1 = subplot(2, 1, 2);
Sp1.Position = [0.2 0.15 0.75 0.32]; % [left bottom width height]
s = shadedErrorBar(1e3*Time(1:end-1),  1 * VarCn(1:end-1), Chi2Error_Cn,...
    'transparent', true, 'patchSaturation',0.2);
set(s.edge,'LineWidth',1,'LineStyle','-', 'Color', Red)
s.mainLine.LineWidth = 2;
s.mainLine.Color = Red;
s.patch.FaceColor = Red;
%errorbar(Time(1:end-1), Corr_var_Eq * VarCn(1:end-1), Chi2Error_Cn, 'o')
hold on
plot(1e3*Time(1:end-1), VarUpAnalytical(1,1:end-1), 'k--', LineWidth = 3)
axis tight
xlim([-2, 11])
box on
ylabel('$\langle x^2(t) \rangle \rm{ ~ [m^2]}$', 'interpreter', 'latex')
xlabel('$t  ~\rm{[ms]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 20, 'fontname', 'times')
set(gcf, 'units', 'centimeters', 'position', [1 1 12 12]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "ExpVar.png",'Resolution', 200)
disp(' ============================================== ')
%% Check that we are within the linear regime
disp('Check: linearity of the perturbation')
VarObservable = mean(var(BetaFitPsd_Eq^2 * EnsembleCn.^2, 0, 2));
figure()
for k = 2:10
    plot(1e3*Time, 1e12* BetaFitPsd_Eq^2 * EnsembleCn(:,k).^2, 'color', [.8 .8 .8], 'LineWidth', 1)
    hold on
end
plot(1e3*Time(1:end-1), 1e12 * (VarCn(1:end-1) - Sf_Cn), 'color', Red, 'LineWidth', 3)
yline(1e12*sqrt(VarObservable), 'k--', 'LineWidth', 1)
axis tight
ylabel('$\{x^2_i(t)\},~ \langle x^2(t) \rangle \rm{ ~ [\mu m^2]}$', 'interpreter', 'latex')
xlabel('$t  ~\rm{[ms]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 20, 'fontname', 'times')
set(gcf, 'units', 'centimeters', 'position', [1 1 12 10]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "CheckLinearRegime.png",'Resolution', 200)
%% steady-state properties
disp('Plot steady-state properties')
Long_eta_s = zeros(numel(Eta_s), 1);
Long_X_s = zeros(numel(X_s), 1);
Long_Xeq_s = zeros(numel(X_s), 1);
for k = 1:numel(Eta_s(1,:))
     Long_eta_s((k-1) * numel(Eta_s(:, 1)) + 1 : k*numel(Eta_s(:, 1))) = Eta_s(:, k);
     Long_X_s((k-1) * numel(X_s(:, 1)) + 1 : k*numel(X_s(:, 1))) = X_s(:, k);
     Long_Xeq_s((k-1) * numel(Xeq_s(:, 1)) + 1 : k*numel(Xeq_s(:, 1))) = Xeq_s(:, k);
end
Long_eta_s = Long_eta_s - mean(Long_eta_s);
Long_X_s = Long_X_s - mean(Long_X_s);
Long_Xeq_s = Long_Xeq_s - mean(Long_Xeq_s);
% plot distributions
[H_Cn, Edge_Cn] = histcounts(Long_X_s, 100, 'Normalization', 'pdf');
[H_Eq, Edge_Eq] = histcounts(Long_Xeq_s, 100, 'Normalization', 'pdf');
figure()
plot(1e6*Edge_Eq(1:end-1), 1e-6*H_Eq, 'Color', Blue, 'LineWidth', 3)
hold on
LongEdge = linspace(min(Edge_Cn), -min(Edge_Cn), 1000);
plot(1e6*LongEdge, 1e-6*sqrt(Kf_Equipart / (2 * pi * kB * T)) ...
        .* exp(-Kf_Equipart .* LongEdge.^2 / (2 * kB*T)), '--', 'Color', 'k', 'LineWidth', 2)

plot(1e6*Edge_Cn(1:end-1), 1e-6*H_Cn, 'Color', Red, 'LineWidth', 3)
axis tight
Sf_Cn_th = D / Omega_f + SigmaEta^2 / (Omega_f * gamma^2 * (Omega_f + OmegaC));
plot(1e6*LongEdge, 1e-6*sqrt(1 / (2 * pi * Sf_Cn_th)) ...
        .* exp(-LongEdge.^2 / (2 * Sf_Cn_th)), '--', 'Color', 'k', 'LineWidth', 2)
ylabel('$P(x) ~\rm{[\mu m^{-1}]}$', 'interpreter', 'latex')
xlabel('$x  ~\rm{[\mu m]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 20, 'fontname', 'times')
set(gcf, 'units', 'centimeters', 'position', [1 1 12 8]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "PDFs.png",'Resolution', 200)
% plot 2D histogram
Data_X_Noise = [Long_X_s, Long_eta_s];
NBins = 80;
bin_x = linspace(min(Long_X_s)/2, -min(Long_X_s)/2, NBins);
bin_y = linspace(min(Long_eta_s)/2, -min(Long_eta_s)/2, NBins);
figure()
cnt = hist3(Data_X_Noise,{bin_x,bin_y});
% Plot frequency values with surf
[x, y] = meshgrid(bin_x, bin_y);
imagesc(1e6*bin_x, 1e12*bin_y, cnt')
axis xy
hold on
%view(2)
n = 256;  % total number of colors
FirstPart = floor(1 * n / 4);
SecondPart = n - FirstPart;
r1 = linspace( 1 , .93, FirstPart)';
g1 = linspace(.95, .68,  FirstPart)';
b1 = linspace(.9, .53, FirstPart)';

r2 = linspace(.93, .9, SecondPart)';
g2 = linspace(.68, .3, SecondPart)';
b2 = linspace(.53, .1, SecondPart)';

colorMap = [r1 g1 b1; r2 g2 b2];

colormap(colorMap)
hold on
 
StartPlot = 11800;
Nplot = 400;
xx = 1e6 * smoothdata(Long_X_s(StartPlot:StartPlot+Nplot), 'movmean', 10);
yy = 1e12 * smoothdata(Long_eta_s(StartPlot:StartPlot+Nplot), 'movmean', 10);
% uncomment to also show a trajectory in the noise-position space
%p1 = plot(xx, yy, ...
%    '-', 'color', Gold, 'LineWidth', 1);
% Cmap3 = customcolormap([0, 1], [StrongRed; Yellow], Nplot);
% for i = 1:Nplot
%     p1 = plot(xx(i), yy(i), 'o', 'Color', Cmap3(i,:));
%     p1.MarkerSize = 3;
%     p1.MarkerFaceColor = Cmap3(i,:);
%     hold on
% end
xlim([-.1, .1])
ylim([-2, 2])
xlabel('$ x  ~\rm{[\mu m]}$', 'fontname', 'times', 'Interpreter', 'latex')
ylabel('$\eta ~\rm{[pN]}$', 'fontname', 'times', 'Interpreter', 'latex')
box on
set(gca,'FontSize', 20, 'fontname', 'times');
set(gcf, 'units', 'centimeters', 'position', [1 1 12 8]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "TwoDpdf.png",'Resolution', 200)
% PSD of noise and trajectories
Omega_X = 2 * pi * f_X;
Omega_Eta = 2 * pi * f_eta;
fminFFT = 0.5;
Nfft = floor(facq/fminFFT);
[S_X, f_X] = pwelch(SteadyTraj(1,:),  hanning(Nfft), floor(Nfft*0.8), Nfft, facq, 'twosided');
figure()
loglog(f_X, S_X, 'o', 'Color', Blue, 'LineWidth', 1)
hold on
[S_X, f_X] = pwelch(SteadyTraj(5,:),  hanning(Nfft), floor(Nfft*0.8), Nfft, facq, 'twosided');
loglog(f_X, S_X, '^', 'Color', Red, 'LineWidth', 1)
loglog(f_X, 2*D./(Omega_X.^2 + Omega_f^2), ...
    '--', 'color', 'k', 'LineWidth', 3)
loglog(f_X, 1./(Omega_X.^2 + Omega_f^2) .* (2*D ...
    + 2 * SigmaEta^2 * OmegaC ./ (gamma^2 * (Omega_X.^2 + OmegaC^2))), ...
    '--', 'color', 'k', 'LineWidth', 3)
axis tight
xlim([1, 1300])
ylabel('$S_{xx}(f) \rm{ ~ [m^2/Hz]}$', 'interpreter', 'latex')
xlabel('$f  ~\rm{[Hz]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 20, 'fontname', 'times')
set(gcf, 'units', 'centimeters', 'position', [1 1 12 8]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "PSDs.png",'Resolution', 200)

% plot trajectories
Time = linspace(0, numel(Eta_s(:,1)) * dt, numel(Eta_s(:,1)));
Pink = [0.84, 0.07, 0.63];
Purple = [0.5, 0.2, 0.66];
Sapin= [0.02, 0.32, 0.26];
figure()
Nplot = 10;
Cmap1 = customcolormap([0 0.5 1],...
    [Pink; OffWhiteWarm; Purple], Nplot);
Cmap2 = customcolormap([0 0.5 1],...
    [Orange; OffWhiteWarm; Red], Nplot);
Sp1 = subplot(2, 1, 1);
Sp1.Position = [0.25 0.57 0.72 0.32]; % [left bottom width height]
for k = 2:Nplot+1
    plot(1e3*Time, 1e12*Eta_s(:,k), 'color', Cmap1(k-1,:), 'LineWidth', 2)
    hold on
end
ylabel('$\eta(t) \rm{ ~ [pN]}$', 'interpreter', 'latex')
xticklabels('')
set(gca, 'fontsize', 18, 'fontname', 'times')
axis tight
ylim([-2, 2])
Sp1 = subplot(2, 1, 2);
Sp1.Position = [0.25 0.2 0.72 0.32]; % [left bottom width height]
for k = 2:Nplot+1
    plot(1e3*Time, 1e6*X_s(:,k), 'color', Cmap2(k-1,:), 'LineWidth', 2)
    hold on
end
axis tight
ylim([-.1, .1])
ylabel('$x(t) \rm{ ~ [\mu m]}$', 'interpreter', 'latex')
xlabel('$t  ~\rm{[ms]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 18, 'fontname', 'times')
set(gcf, 'units', 'centimeters', 'position', [1 1 12 10]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "ExpTraj.png",'Resolution', 200)
disp(' ============================================== ')
%% Compute correlations
disp('Compute correlation function in steady-state')
Vx_s = mean(X_s.^2, 2);
NTime_x = numel(X_s(:,1));
NTraj_x = numel(X_s(1, :));

MaxT_0 =  floor(7*NTime_x/8);
Cxx_eq = Cross_TA_Cxx_flip(Xeq_s, Xeq_s, MaxT_0, NTime_x, NTraj_x);
Cxx = Cross_TA_Cxx_flip(X_s, X_s, MaxT_0, NTime_x, NTraj_x);
Cx2x2 = Cross_TA_Cxx_flip(X_s.^2, X_s.^2, MaxT_0, NTime_x, NTraj_x);
CxEtax2 = Cross_TA_Cxx_flip(X_s .* Eta_s, X_s.^2, MaxT_0, NTime_x, NTraj_x);
CEta2x2 = Cross_TA_Cxx_flip(Eta_s.^2, X_s.^2, MaxT_0, NTime_x, NTraj_x);
% Analytical result
ErrorGamma = 0.05;
LagArray = linspace(dt, numel(Cx2x2)*dt, numel(Cx2x2));
ErrOmega0 = sqrt((ErrKi_PSD/Ki_PSD * Omega_f)^2 + (ErrorGamma * Omega_f)^2);
AnalyticalCxx_Eq = D / Omega_f * exp(-Omega_f * LagArray);

AnalyticalCxx_Cn = D/Omega_f * exp(-Omega_f * LagArray)...
    + SigmaEta_f^2 * OmegaC/ (Omega_f * (OmegaC^2 - Omega_f^2))...
    * ( exp(-Omega_f * LagArray)...
    - Omega_f / OmegaC * exp(-OmegaC * LagArray) );
%% Test one dimensional FRR
disp('Test equilibrium FRR in the 1D position space')
Ntime_i = 350; % until after relax
Ntime_f = Ntime - Ntime_i+1;
VxEq = VarEq(t:Ntime+t-Ntime_i);
Vx = VarCn(t:Ntime+t-Ntime_i);
Time = linspace(0, (Ntime-Ntime_i) * dt, Ntime-Ntime_i+1);
Ki = 13.5e-6;
Kf = 15e-6;
Start = 1;
Stop = min([numel(Cx2x2), numel(VarEq)]);
figure()
Sp1 = subplot(2, 1, 1);
Sp1.Position = [0.2 0.6 0.72 0.32]; % [left bottom width height]
FdrFactor = Kf_Equipart * (Kf_Equipart - Ki_Equipart)/(kB * T * Ki_Equipart); 
plot(1e3*Time(Start:Stop), 1e12 * FdrFactor .* Cxx_eq(Start:Stop).^2, '--', 'color', StrongBlue, 'LineWidth', 4)
hold on
s = shadedErrorBar(Time(Start:Stop) * 1e3, 1e12 * FdrFactor .* Cxx_eq(Start:Stop).^2, ...
    1e12 * 2 * ErrorBeta / BetaFitPsd_Eq * FdrFactor .* Cxx_eq(Start:Stop).^2,...
    'transparent', true, 'patchSaturation',0.2);
set(s.mainLine,'LineWidth', 4,'LineStyle','--', 'Color', StrongBlue)
set(s.edge,'LineWidth', 1,'LineStyle','--', 'Color', Blue)
s.patch.FaceColor = LightBlue;
plot(1e3*Time(Start:Stop), 1e12 * (VxEq(Start:Stop) - Sf_eq), 'color', Blue, 'LineWidth', 3)
s = shadedErrorBar(Time(Start:Stop) * 1e3, 1e12 * (VxEq(Start:Stop) - Sf_eq), ...
    1e12 .* Chi2Error_Eq(Start:Stop),...
    'transparent', true, 'patchSaturation',0.2);
set(s.edge,'LineWidth', 1,'LineStyle','-', 'Color', 'none')
s.mainLine.LineWidth = 2;
s.mainLine.Color = Blue;
s.patch.FaceColor = Blue;
yline(0, 'k-')
axis tight
%ylabel('$\rm{ ~ [\mu m^2]}$', 'interpreter', 'latex')
xticklabels('')
set(gca, 'fontsize', 20, 'fontname', 'times')
Sp1 = subplot(2, 1, 2);
Sp1.Position = [0.2 0.17 0.72 0.34]; % [left bottom width height]
FdrFactor = Kf_Equipart * (Kf_Equipart - Ki_Equipart)/(kB * T * Ki_Equipart); 
plot(1e3*Time(Start:Stop), 1e12 * FdrFactor .* Cxx(Start:Stop).^2, '--', 'color', BrownRed, 'LineWidth', 3)
hold on
s = shadedErrorBar(Time(Start:Stop) * 1e3, 1e12 * FdrFactor .* Cxx(Start:Stop).^2, ...
    1e12 * 2 * ErrorBeta / BetaFitPsd_Eq * FdrFactor .* Cxx(Start:Stop).^2,...
    'transparent', true, 'patchSaturation', 0.2);
set(s.mainLine,'LineWidth', 3,'LineStyle','--', 'Color', BrownRed)
set(s.edge,'LineWidth', 1,'LineStyle','--', 'Color', Red)
s.patch.FaceColor = Orange;
plot(1e3*Time(Start:Stop), 1e12 * (Vx(Start:Stop) - Sf_Cn), 'color', Red, 'LineWidth', 3)
s = shadedErrorBar(Time(Start:Stop) * 1e3, 1e12 * (Vx(Start:Stop) - Sf_Cn), ...
    1e12 .* Chi2Error_Cn(Start:Stop),...
    'transparent', true, 'patchSaturation',0.2);
set(s.edge,'LineWidth', 1,'LineStyle','-', 'Color', Red)
s.mainLine.LineWidth = 2;
s.mainLine.Color = Red;
s.patch.FaceColor = Red;
yline(0, 'k-')
axis tight
xlabel('$t  ~\rm{[ms]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 20, 'fontname', 'times')
set(gcf, 'units', 'centimeters', 'position', [1 1 12 12]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "OneD_Frr.png",'Resolution', 200)
%% Compute the terms in 2D Marokovian-embedding FRR
disp('Compute FRR in 2D Markovian-embedding')
Tau1 = gamma/Kf_Equipart;
tt = Tau1 + TauC;
M = D * gamma^2 * tt^2 + D2 * Tau1^2 * TauC^2;
a0 = (D * gamma^2 * tt^3 + D2 * Tau1^2 * TauC^2 * (Tau1 + 3 * TauC)) ...
        / (2 * tt * M );
a11 = - (gamma^2 * tt * (D * gamma^2 * tt^3 + D2 * Tau1^2 * TauC^2 * (Tau1 + 3 * TauC))) ...
        / (2 *Tau1 * M^2);
a12 = (gamma * Tau1 * TauC * (D * gamma^2 * tt^2 + D2 * Tau1 * TauC^2 * (Tau1 + 2 * TauC))) ...
        / (M^2);
a22 = - (D * gamma^2 * Tau1 * TauC^2 * (Tau1^2 - TauC^2) + D2 * Tau1^3 * TauC^4) ...
        / (2* M^2);

save('Time_exp.mat', "Time");
save('Vx_exp.mat', "Vx");
save('Cx2x2_exp.mat', "Cx2x2");
save('CxEtax2_exp.mat', "CxEtax2");
save('CEta2x2_exp.mat', "CEta2x2");

Start = 1;
Stop = min([numel(Cx2x2), numel(Vx)]);
figure()
subplot(2,2,1)
plot(Time(Start:Stop), Vx(Start:Stop), 'LineWidth', 2);
ylabel('$\langle x^2 \rangle \rm{[m^2]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 15, 'fontname', 'times')
axis tight
subplot(2,2,2)
plot(Time(Start:Stop), a11 * Cx2x2(Start:Stop), 'LineWidth', 2);
ylabel('$a_{11} \langle x^2(0) x^2(t) \rangle_s \rm{[m^2]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 15, 'fontname', 'times')
axis tight
subplot(2,2,3)
plot(Time(Start:Stop), a12 * CxEtax2(Start:Stop), 'LineWidth', 2);
ylabel('$a_{12} \langle x(0) \eta(0) x^2(t) \rangle_s \rm{[m^4]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 15, 'fontname', 'times')
axis tight
subplot(2,2,4)
plot(Time(Start:Stop), a22 * CEta2x2(Start:Stop), 'LineWidth', 2);
ylabel('$a_{22} \langle \eta^2(0) x^2(t) \rangle_s \rm{[m^2]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 15, 'fontname', 'times')
axis tight
set(gcf, 'units', 'centimeters', 'position', [1 1 25 18]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "ExpAllCorr.png",'Resolution', 200)
%% Plot two-dim FDR
Y = (a0 * Vx_s(1:numel(Cx2x2))' + a11 * Cx2x2 + a12 * CxEtax2 + a22 * CEta2x2);
FRR_Rhs = (Ki_Equipart - Kf_Equipart) / Kf_Equipart * Y;
Err_FRR_Rhs = sqrt( (ErrKi_PSD/Kf)^2 + (ErrKf_PSD/Kf)^2 + ((Ki - Kf)/Kf^2 * ErrKf_PSD)^2 ) * Y;
FRR_Lhs = Vx - Sf_Cn;

save('FRR_Lhs_exp.mat', "FRR_Lhs");
save('FRR_Rhs_exp.mat', "FRR_Rhs");
save('Err_FRR_Rhs_exp.mat', "Err_FRR_Rhs");

figure()
Start = 2;
semilogx(Time(Start:end)*1e3, 1e12 * FRR_Lhs(Start:end), '-', 'LineWidth', 3, 'color', Red)
hold on
s = shadedErrorBar(Time(Start:end) * 1e3, 1e12 * (Vx(Start:end) - Sf_Cn), 1e12 * Chi2Error_Cn(Start:numel(Vx)),...
    'transparent', true, 'patchSaturation',0.2);
set(s.edge,'LineWidth', 1,'LineStyle','-', 'Color', Red)
s.mainLine.LineWidth = 2;
s.mainLine.Color = Red;
s.patch.FaceColor = Red;
hold on
plot(Time(Start:numel(FRR_Rhs))*1e3, 1e12 * FRR_Rhs(Start:end), '--', 'LineWidth', 4, 'Color', BrownRed)
s = shadedErrorBar(Time(Start:numel(FRR_Rhs)) * 1e3, 1e12 * FRR_Rhs(Start:end), 1e12 * 2 * ErrorBeta / BetaFitPsd_Eq * FRR_Rhs(Start:end),...
    'transparent', true, 'patchSaturation',0.2);
set(s.mainLine,'LineWidth', 4,'LineStyle','--', 'Color', BrownRed)
set(s.edge,'LineWidth', 2,'LineStyle','--', 'Color', BrownRed)
s.patch.FaceColor = Orange;
yline(0, 'k-', 'LineWidth', 1)
axis tight
xlim([0, max(Time(Start:numel(FRR_Rhs))*1e3)])
xticks([.1, 1, 10])
xticklabels({'0.1', '1', '10'})
box on
ylabel('$\rm{FRR ~ [\mu m^2]}$', 'interpreter', 'latex')
xlabel('$t  ~\rm{[ms]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 20, 'fontname', 'times')
set(gcf, 'units', 'centimeters', 'position', [1 1 12 12]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "ExpFRR.png",'Resolution', 200)
disp(' ============================================== ')
%% clear ensembles
clear EnsembleEq;
clear EnsembleCn;
clear EnsembleNoiseCn;
%% Numerlical simlulations
disp('Numerical simulation of Langevin equation uing the experimental parameters')
Tau1_i = gamma / Ki_Equipart;
Tau1_f = gamma / Kf_Equipart;
TotalTime_sim = 12 * Tau1_f;
Tau2 = 1./OmegaC;
dt_sim = 1e-5;
NTime_sim = floor(TotalTime_sim / dt_sim);
NTraj_sim = 150000;
% Generate trajecotories in steady-state
InitialCondition = sqrt(SigmaEta^2) .* random('Normal', 0, 1, [NTraj_sim, 1]);
Eta_s_sim = OrnUhlFunction_noisyStiffness(NTraj_sim, NTime_sim, dt_sim, D2, 1, 1/Tau2 * ones([NTraj_sim, NTime_sim]), zeros([NTraj_sim, NTime_sim]), InitialCondition);
InitialCondition = sqrt(D * Tau1_i + D2 * Tau1_i^2 * Tau2^2 / (gamma^2 * (Tau1_i + Tau2))) .* random('Normal', 0, 1, [NTraj_sim, 1]);
X_s_sim = OrnUhlFunction_noisyStiffness(NTraj_sim, NTime_sim, dt_sim, D, gamma, Kf_Equipart * ones([NTraj_sim, NTime_sim]), Eta_s_sim, InitialCondition);

% Generate trajecotories with stiffness change
NTime_i_sim = floor(NTime_sim / 3);
Kappa = cat(2, Ki_Equipart * ones([NTraj_sim, NTime_i_sim]), Kf_Equipart * ones([NTraj_sim, NTime_sim]));
InitialCondition = sqrt(SigmaEta^2) .* random('Normal', 0, 1, [NTraj_sim, 1]);
Eta_sim = OrnUhlFunction_noisyStiffness(NTraj_sim, NTime_sim, dt_sim, D2, 1, 1/Tau2 * ones([NTraj_sim, NTime_sim]), zeros([NTraj_sim, NTime_sim]), InitialCondition);
InitialCondition = sqrt(D * Tau1_i + D2 * Tau1_i^2 * Tau2^2 / (gamma^2 * (Tau1_i + Tau2))) .* random('Normal', 0, 1, [NTraj_sim, 1]);
X_sim = OrnUhlFunction_noisyStiffness(NTraj_sim, NTime_sim, dt_sim, D, gamma, Kappa, Eta_sim, InitialCondition);
X_sim(:, 1:NTime_i_sim) = []; % cut the beginning
Eta_sim(:, 1:NTime_i_sim) = [];
X_s_sim(:, 1:NTime_i_sim) = []; % cut the beginning
Eta_s_sim(:, 1:NTime_i_sim) = [];
NTime_x_sim = numel(X_sim(1,:));
%% Compute variances and correlations
disp('Compute correlation of the trajectories resulting from numerical simulation')
Vx_sim = var(X_sim, 0, 1);
MaxT_0_sim =  floor(3*NTime_x_sim/5);
NTraj_Corr_sim = floor(NTraj_sim / 50); % smaller statistics is enought for time-ensemble average
Cxx_sim = Cross_TA_Cxx(X_s_sim.^2, X_s_sim.^2, MaxT_0_sim, NTime_x_sim, NTraj_Corr_sim);
CxEtaX_sim = Cross_TA_Cxx(X_s_sim .* Eta_s_sim, X_s_sim.^2, MaxT_0_sim, NTime_x_sim, NTraj_Corr_sim);
CxEta_sim = Cross_TA_Cxx(Eta_s_sim.^2, X_s_sim.^2, MaxT_0_sim, NTime_x_sim, NTraj_Corr_sim);
Time_sim = linspace(0, NTime_x_sim * dt_sim, NTime_x_sim);
%% plot only cross corr with experiment
Start = 1;
Stop_sim = numel(Cxx_sim);
figure()
subplot(2,2,1)
plot(1e3*Time_sim(Start:Stop_sim), Vx_sim(Start:Stop_sim) - Vx_sim(Stop_sim), '-.', 'Color', Purple, 'LineWidth', 3);
hold on
plot(1e3*Time(1:numel(Vx)), Vx - Sf_Cn, 'Color', Red, 'LineWidth', 3);
%ylabel('$a_{11} \langle x^2(0) x^2(t) \rangle_s \rm{[m^2]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 18, 'fontname', 'times')
axis tight
xticklabels('')

subplot(2,2,2)
plot(1e3*Time_sim(Start:Stop_sim), a11 * Cxx_sim(Start:Stop_sim), '-.', 'Color', Purple, 'LineWidth', 3);
hold on
plot(1e3*Time(1:numel(Cx2x2)), a11 * Cx2x2, 'Color', Red, 'LineWidth', 3);
%ylabel('$a_{11} \langle x^2(0) x^2(t) \rangle_s \rm{[m^2]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 18, 'fontname', 'times')
axis tight
xticklabels('')
subplot(2,2,3)
plot(1e3*Time_sim(Start:Stop_sim), a12 * CxEtaX_sim(Start:Stop_sim), '-.', 'Color', Purple, 'LineWidth', 3);
hold on 
plot(1e3*Time(1:numel(CxEtax2)), a12 * CxEtax2, 'Color', Red, 'LineWidth', 3);
%ylabel('$a_{12} \langle x(0) \eta(0) x^2(t) \rangle_s \rm{[m^4]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 18, 'fontname', 'times')
axis tight
xlabel('$t~ \rm{[ms]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 18, 'fontname', 'times')
subplot(2,2,4)
plot(1e3*Time_sim(Start:Stop_sim), a22 * CxEta_sim(Start:Stop_sim), '-.', 'Color', Purple, 'LineWidth', 3);
hold on 
plot(1e3*Time(1:numel(CEta2x2)), a22 * CEta2x2, 'Color', Red, 'LineWidth', 3);
%ylabel('$a_{22} \langle \eta^2(0) x^2(t) \rangle_s \rm{[m^2]}$', 'interpreter', 'latex')
xlabel('$t~ \rm{[ms]}$', 'interpreter', 'latex')
set(gca, 'fontsize', 18, 'fontname', 'times')
axis tight
set(gcf, 'units', 'centimeters', 'position', [1 1 25 12]);
set(gcf, 'Color', 'w');
exportgraphics(gcf, "SimCrossCorr.png",'Resolution', 200)

