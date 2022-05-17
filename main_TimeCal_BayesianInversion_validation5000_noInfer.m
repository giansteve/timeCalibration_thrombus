% Bayesian inversion of the Time Calibration project

% Setting the environment
clearvars
close all
clc

% Graphic settings
set(0,'DefaultFigureWindowStyle','default')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex')
set(0,'defaultAxesFontSize',11)

% ROOT Destination
cd C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\timeCalibration_thrombus
root_destination = pwd;
addpath(root_destination)

% initialize UQlab
% uqlab

%% Load RawData
load('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_valid5000noInfer_final\RawOutput_04May22_5000sim.mat')

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

%% plot fitting MRI
folderPath = 'TimeCal_MRIdata';
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel_Validation5000_noInfer\\%s',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotFittingMRI(fittedData,human_thr,'MRI_fitting')
cd(root_destination)


%% Load output: H/S and L/S
phic_HS_threshold = outToPCE.H_S;
phic_LS_threshold = outToPCE.L_S;

%% Define Input variables names
var_names = {'${D}_{\mathrm{c}}$',...
    '$\langle{\dot{\gamma}}\langle_\mathrm{t} $'};

%% Statistics on the input
fprintf('Get statistics INPUT ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel_Validation5000_noInfer\\StatsInput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
figure('Visible','off')
plotmatrix(exp_design(:,[1 6]))
GM_printBMP(600,600,'InputStats')
GM_printEPS(600,600,'InputStats')
close all
cd(root_destination)

%% Statistics on the output
fprintf('Get statistics OUTPUT ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel_Validation5000_noInfer\\StatsOutput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_statisticsOutput(phic_HS_threshold,phic_LS_threshold,human_thr);
cd(root_destination)





















