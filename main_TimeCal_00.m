% Main script for Time Calibration of the thrombus model. Adjustment
% started on 08/12/21 with the scope of making the code faster.

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
cd M:\IFM\User\melito\Server\Projects\TimeCalibration
root_destination = pwd;
addpath(root_destination)

% initialize UQlab
% uqlab

% consider adding them in the folder
% addpath('M:\IFM\User\melito\Server\Projects\matlab_funct\SA')
% addpath('M:\IFM\User\melito\Server\Projects\matlab_funct\general')
% addpath('M:\IFM\User\melito\Server\Projects\matlab_funct')

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
    dest_plot = sprintf('Plot\\%s',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotFittingMRI(fittedData,human_thr,'MRI_fitting')
cd(root_destination)

%% Define Input variables names
var_names = {'${D}_{\mathrm{c}}$','${k}_{\mathrm{c}}$','${k}_{\mathrm{{BP}}}$','${c}_{\mathrm{t}}$',...
    '${c}_{\mathrm{{BPt}}}$','$\overline{T}_{\mathrm{Rt}}$','${k}_{\mathrm{{c,wall}}}$',...
    '${c}_{\mathrm{{AP}}}$','${c}_{\mathrm{{BPbt}}}$','$\dot{\overline{\gamma}}_t$'};

%% Read input file
myVars = {'INPUT','exp_design'};
load('timeCal_00input_2000sims_4Feb21.mat',myVars{:})

%% Load output: H/S and L/S
load('TimeCal_OUTPUT.mat')
% delete crushed sim from the input
exp_design(crushed_sim_idx,:) = [];
Ns = size(exp_design,1);
M = size(exp_design,2);

%% Statistics on the input
fprintf('Get statistics INPUT ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot\\StatsInput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
figure('Visible','off')
plotmatrix(exp_design)
GM_printBMP(600,600,'InputStats')
GM_printEPS(600,600,'InputStats')
close all
cd(root_destination)

%% Statistics on the output
fprintf('Get statistics OUTPUT ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot\\StatsOutput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_statisticsOutput(OUT_HS,OUT_LS,OUT_time,fittedData);
cd(root_destination)

%% Create metamodel of the thrombus model
% moving PCE computation to UQ[py]lab
MRI_time_index = (1:10:61);
metamodel.Type = 'Metamodel';
metamodel.MetaType = 'PCE';
metamodel.Display = 'verbose';
metamodel.Method = 'lars';
metamodel.Degree = 3:12;
metamodel.TruncOptions.qNorm = 0.95;
metamodel.DegreeEarlyStop = false;
metamodel.Input = INPUT;
metamodel.ExpDesign.NSamples = Ns;
metamodel.ExpDesign.X = exp_design;
% PCE H/S
metamodel.ExpDesign.Y = OUT_HS(MRI_time_index,:)';
PCE_HS = uq_createModel(metamodel);
% PCE L/S
metamodel.ExpDesign.Y = OUT_LS(MRI_time_index,:)';
PCE_LS = uq_createModel(metamodel);

% save PS: remember to create function to save in specific folder with the
% date :)
save('TimeCal_postSurrogate_101221')

%% Residuals MRI (fitted data) and MODEL
% static residuals: only at the data points
% static_time_index = (1:10:61);
% ns = size(human_thr.H_S,1); % size of the residual points
% r_s_HS = (1/sqrt(ns)).*( OUT_HS(static_time_index,:) - human_thr.H_S(:,1) )./human_thr.H_S(:,1);
% r_s_HS(1,:) = 0; % imposing first time instance with 0 discrepancy = 0/0 -> NaN
% % dynamic residuals: for the whole process of thrombus formation
% r_d_HS = (1/sqrt(size(fittedData.H_S,2))).*( OUT_HS - fittedData.H_S' )./fittedData.H_S';
% r_d_HS(1,:) = 0; % imposing first time instance with 0 discrepancy = 0/0 -> NaN













































