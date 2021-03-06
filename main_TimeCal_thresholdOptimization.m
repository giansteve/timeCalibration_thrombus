% TimeCal phi_c optimization algorithm

% Setting the environment
% clearvars
% close all
% clc

% Graphic settings
set(0,'DefaultFigureWindowStyle','default')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex')
set(0,'defaultAxesFontSize',11)

% ROOT Destination
cd M:\IFM\User\melito\Server\Projects\TimeCalibration
root_destination = pwd;
addpath(root_destination)

% initialize UQlab
% uqlab

% load raw data from model sim
% load('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\TimeCal_phic.mat')
% load MRI data and model sim rest of data
% load('TimeCal_OUTPUT.mat')

%% Dealing with Experimental Data
% ===============================SOURCE===================================
% Taylor, Joshua O. et al. "Development of a computational model for
% macroscopic predictions of device-induced thrombosis", J. Biomechanics
% and Modeling in Mechanobiology, 2016
% ========================================================================
% Human thrombus in real-time
human_thr.time = xlsread('MRIdata.xlsx',1,'B20:B26');                       % time
human_thr.volume = xlsread('MRIdata.xlsx',1,'C20:C26');                     % thrombus volume
human_thr.volume(:,2) = xlsread('MRIdata.xlsx',1,'D20:D26');                % std thrombus volume
human_thr.exposedArea = xlsread('MRIdata.xlsx',1,'E20:E26');                % surface area
human_thr.exposedArea(:,2) = xlsread('MRIdata.xlsx',1,'F20:F26');           % std surface area
human_thr.H_S = xlsread('MRIdata.xlsx',1,'G20:G26');                        % Height/Step
human_thr.H_S(:,2) = xlsread('MRIdata.xlsx',1,'H20:H26');                   % std Height/Step
human_thr.L_S = xlsread('MRIdata.xlsx',1,'I20:I26');                        % Length/Step
human_thr.L_S(:,2) = xlsread('MRIdata.xlsx',1,'J20:J26');                   % std Length/Step

%% Fitting data into 'exp2' function (educated guess)
% function_type = a*exp(b*x) + c*exp(d*x) = exp2
[fitData.ff_HS,fitData.gof_HS] = fit(human_thr.time,human_thr.H_S(:,1),'exp2');
[fitData.ff_LS,fitData.gof_LS] = fit(human_thr.time,human_thr.L_S(:,1),'exp2');
[fitData.ff_SA,fitData.gof_SA] = fit(human_thr.time,human_thr.exposedArea(:,1),'exp2');
[fitData.ff_Vol,fitData.gof_Vol] = fit(human_thr.time,human_thr.volume(:,1),'exp2');

%% Generate new data points based on previous fit
n_pt = 61;                                                                      % nr. pts to generate
fittedData.time = linspace(min(human_thr.time),max(human_thr.time),n_pt);       % time [min]
fittedData.H_S = fitData.ff_HS.a*exp(fitData.ff_HS.b.*fittedData.time) + ...    % Height/Step
    fitData.ff_HS.c*exp(fitData.ff_HS.d.*fittedData.time);
fittedData.L_S = fitData.ff_LS.a*exp(fitData.ff_LS.b.*fittedData.time) + ...    % Length/Step
    fitData.ff_LS.c*exp(fitData.ff_LS.d.*fittedData.time);
fittedData.SA = fitData.ff_SA.a*exp(fitData.ff_SA.b.*fittedData.time) + ...     % surface area
    fitData.ff_SA.c*exp(fitData.ff_SA.d.*fittedData.time);
fittedData.Vol = fitData.ff_Vol.a*exp(fitData.ff_Vol.b.*fittedData.time) + ...  % thrombus volume
    fitData.ff_Vol.c*exp(fitData.ff_Vol.d.*fittedData.time);
% delete negative and incoherent data
fittedData.H_S(fittedData.H_S < 0) = 0;
fittedData.L_S(fittedData.L_S < 0) = 0;
fittedData.SA(fittedData.SA < 0) = 0;
fittedData.Vol(fittedData.Vol < 0) = 0;
% transformation into SI (from cm to m)
human_thr.exposedArea = human_thr.exposedArea ./ 10000;
human_thr.volume = human_thr.volume ./ 1000000;
fittedData.SA = fittedData.SA ./ 10000;
fittedData.Vol = fittedData.Vol ./ 1000000;
% post-processing: get growth rate of each quantity
deltaT = fittedData.time(2);
human_thr.growthRate.H_S = diff(fittedData.H_S)/(deltaT);           % 1/min
human_thr.growthRate.L_S = diff(fittedData.L_S)/(deltaT);           % 1/min
human_thr.growthRate.volume = diff(fittedData.Vol)/(deltaT);        % m^3/min
human_thr.growthRate.exposedArea = diff(fittedData.SA)/(deltaT);    % m^2/min

%% Optimization starts
threshold_val = 0.005:0.005:0.2;

% new fitted data
n_pt = 3000;                                                                      % nr. pts to generate
MRIfit_optim.time = linspace(min(human_thr.time),max(human_thr.time),n_pt);       % time [min]
MRIfit_optim.H_S = fitData.ff_HS.a*exp(fitData.ff_HS.b.*MRIfit_optim.time) + ...    % Height/Step
    fitData.ff_HS.c*exp(fitData.ff_HS.d.*MRIfit_optim.time);
MRIfit_optim.L_S = fitData.ff_LS.a*exp(fitData.ff_LS.b.*MRIfit_optim.time) + ...    % Length/Step
    fitData.ff_LS.c*exp(fitData.ff_LS.d.*MRIfit_optim.time);

% H/S, totPoints: 50, lenght(2):0.0025
% L/S, totPoints: 160, lenght(2):0.0695
n_sim_toOptimize = size(OUT_HS,2); % number of succeded simulations
% pre-allocation of memory
HS_toOptim = zeros(3000,size(threshold_val,2));
LS_toOptim = zeros(3000,size(threshold_val,2));
thresh_opt_val_HS = zeros(1,3000);
thresh_opt_avg_HS = zeros(n_sim_toOptimize,1);
thresh_opt_val_LS = zeros(1,3000);
thresh_opt_avg_LS = zeros(n_sim_toOptimize,1);
% begin optimization of thrombus threshold value
for nsim = 1:n_sim_toOptimize
    for tresh_idx = 1:size(threshold_val,2)
        logic_temp = phic_HS(:,:,nsim) >= threshold_val(tresh_idx); % count cell thrombus above threshold H/S
        HS_toOptim(:,tresh_idx) = (sum(logic_temp,2) / 50 * 0.0025) / 0.0025; % get corresponding H/S
        logic_temp = phic_LS(:,:,nsim) >= threshold_val(tresh_idx); % count cell thrombus above threshold L/S
        LS_toOptim(:,tresh_idx) = (sum(logic_temp,2) / 160 * 0.0695) / 0.0025; % get corresponding L/S
    end
    % determine residual: Model-MRI
    % H/S
    res_HS = ((MRIfit_optim.H_S' - HS_toOptim).^2);
    [~,midx2] = min(res_HS,[],2);
    thresh_opt_val_HS = threshold_val(midx2); % optimum threshold in time
    thresh_opt_avg_HS(nsim) = mean(thresh_opt_val_HS); % time avg optimum threshold
    % L/S
    res_LS = ((MRIfit_optim.L_S' - LS_toOptim).^2);
    [~,midx2] = min(res_LS,[],2);
    thresh_opt_val_LS = threshold_val(midx2); % optimum threshold in time
    thresh_opt_avg_LS(nsim) = mean(thresh_opt_val_LS); % time avg optimum threshold
    
    %     % plot for verification (only H/S)
    %     figure(1)
    %     subplot(211)
    %     plot(linspace(0,1,3000),LS_toOptim)
    %     hold on
    %     plot(linspace(0,1,7),human_thr.L_S(:,1),'k.','MarkerSize',7)
    %     plot(linspace(0,1,3000),MRIfit_optim.L_S,'k-')
    %     ylabel('H/S')
    %     hold off
    %
    %     subplot(212)
    %     plot(linspace(0,1,3000),thresh_opt_val_LS,'r')
    %     hold on
    %     yline(thresh_opt_avg_LS(nsim),'k');
    %     set(gca,'YScale','linear')
    %     hold off
    %     xlabel('time')
    %     ylabel('threshold val')
    
end

%% Plot results of optimum threshold
threshold_edges = threshold_val;
centres = threshold_edges(1:end-1)+ diff(threshold_edges)/2;

cd(root_destination)
try
    dest_plot = sprintf('Plot\\StatsOutput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
figure('Visible','off')
subplot(211)
h_HS = histcounts(thresh_opt_avg_HS,threshold_edges,'Normalization','probability');
plot(centres,(smooth(h_HS)),'k')
title('H/S')
grid on
ylabel('prob minimum residual')
subplot(212)
h_LS = histcounts(thresh_opt_avg_LS,threshold_edges,'Normalization','probability');
plot(centres,(smooth(h_LS)),'k')
title('L/S')
grid on
xlabel('threshold val')
ylabel('prob minimum residual')
GM_printBMP(400,400,'threshold_optimization')
GM_printEPS(400,400,'threshold_optimization')
cd(root_destination)

%% Select threshold value
thrombus_threshold = 0.01; % threshold selected
% pre-allocation
phic_HS_threshold = zeros(3000,n_sim_toOptimize);
phic_LS_threshold = phic_HS_threshold;
for nsim = 1:n_sim_toOptimize
    logic_temp = phic_HS(:,:,nsim) >= thrombus_threshold; % count cell thrombus above threshold H/S
    phic_HS_threshold(:,nsim) = (sum(logic_temp,2) / 50 * 0.0025) / 0.0025; % get corresponding H/S
    logic_temp = phic_LS(:,:,nsim) >= thrombus_threshold; % count cell thrombus above threshold L/S
    phic_LS_threshold(:,nsim) = (sum(logic_temp,2) / 160 * 0.0695) / 0.0025; % get corresponding L/S
end

cd(root_destination)
try
    dest_plot = sprintf('Plot\\StatsOutput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
figure('Visible','off')
subplot(211)
plot(linspace(0,1,3000),phic_HS_threshold,'Color',[.5 .5 .5])
hold on
plot(linspace(0,1,7),human_thr.H_S(:,1),'k.','MarkerSize',7)
plot(linspace(0,1,3000),MRIfit_optim.H_S,'k-')
ylabel('$H/S ~(\%)$')
subplot(212)
plot(linspace(0,1,3000),phic_LS_threshold,'Color',[.5 .5 .5])
hold on
plot(linspace(0,1,7),human_thr.L_S(:,1),'k.','MarkerSize',7)
plot(linspace(0,1,3000),MRIfit_optim.L_S,'k-')
ylabel('$L/S ~(\%)$')
xlabel('$t^*$')
GM_printBMP(400,400,'threshold_application')
GM_printEPS(400,400,'threshold_application')
cd(root_destination)

% save outputs
save('TimeCal_OUTPUT_allTime.mat','phic_HS_threshold','phic_LS_threshold')



