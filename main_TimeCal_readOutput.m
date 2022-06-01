% use this file to read the output generated from the model simulations

clearvars
close all
clc

set(0,'DefaultFigureWindowStyle','default')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex')
set(0,'defaultAxesFontSize',11)

% cd M:\IFM\User\melito\Server\Projects\TimeCalibration
root_destination = pwd;
addpath(root_destination)

addpath('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\matlab_funct\SA')
addpath('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\matlab_funct\general')
addpath('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\matlab_funct')

%% Experimental Data
% Human real-time
% human_thr.time = xlsread('MRIdata.xlsx',1,'B20:B26');                       % time
% human_thr.volume = xlsread('MRIdata.xlsx',1,'C20:C26');                     % thrombus volume
% human_thr.volume(:,2) = xlsread('MRIdata.xlsx',1,'D20:D26');                % std thrombus volume
% human_thr.exposedArea = xlsread('MRIdata.xlsx',1,'E20:E26');                % thrombus area
% human_thr.exposedArea(:,2) = xlsread('MRIdata.xlsx',1,'F20:F26');           % std thrombus area
% human_thr.H_S = xlsread('MRIdata.xlsx',1,'G20:G26');                        % Height/Step
% human_thr.H_S(:,2) = xlsread('MRIdata.xlsx',1,'H20:H26');                   % std Height/Step
% human_thr.L_S = xlsread('MRIdata.xlsx',1,'I20:I26');                        % Length/Step
% human_thr.L_S(:,2) = xlsread('MRIdata.xlsx',1,'J20:J26');                   % std Length/Step

%% Function fitting
% [fitData.ff_HS,fitData.gof_HS] = fit(human_thr.time,human_thr.H_S(:,1),'exp2');
% [fitData.ff_LS,fitData.gof_LS] = fit(human_thr.time,human_thr.L_S(:,1),'exp2');
% [fitData.ff_SA,fitData.gof_SA] = fit(human_thr.time,human_thr.exposedArea(:,1),'exp2');
% [fitData.ff_Vol,fitData.gof_Vol] = fit(human_thr.time,human_thr.volume(:,1),'exp2');
% % fitted data
% % function_type = a*exp(b*x) + c*exp(d*x);
% fittedData.time = linspace(min(human_thr.time),max(human_thr.time),61);    % time [min]
% fittedData.H_S = fitData.ff_HS.a*exp(fitData.ff_HS.b.*fittedData.time) + fitData.ff_HS.c*exp(fitData.ff_HS.d.*fittedData.time);
% fittedData.L_S = fitData.ff_LS.a*exp(fitData.ff_LS.b.*fittedData.time) + fitData.ff_LS.c*exp(fitData.ff_LS.d.*fittedData.time);
% fittedData.SA = fitData.ff_SA.a*exp(fitData.ff_SA.b.*fittedData.time) + fitData.ff_SA.c*exp(fitData.ff_SA.d.*fittedData.time);
% fittedData.Vol = fitData.ff_Vol.a*exp(fitData.ff_Vol.b.*fittedData.time) + fitData.ff_Vol.c*exp(fitData.ff_Vol.d.*fittedData.time);
% fittedData.H_S(fittedData.H_S < 0) = 0;
% fittedData.L_S(fittedData.L_S < 0) = 0;
% fittedData.SA(fittedData.SA < 0) = 0;
% fittedData.Vol(fittedData.Vol < 0) = 0;
% % transformation into SI
% human_thr.exposedArea = human_thr.exposedArea ./ 10000;
% human_thr.volume = human_thr.volume ./ 1000000;
% fittedData.SA = fittedData.SA ./ 10000;
% fittedData.Vol = fittedData.Vol ./ 1000000;
% % post-processing experimental fitted data
% human_thr.growthRate.H_S = diff(fittedData.H_S)/(fittedData.time(2));       % 1/min
% human_thr.growthRate.L_S = diff(fittedData.L_S)/(fittedData.time(2));       % 1/min
% human_thr.growthRate.volume = diff(fittedData.Vol)/(fittedData.time(2));    % m^3/min
% human_thr.growthRate.exposedArea = diff(fittedData.SA)/(fittedData.time(2));% m^2/min

%% Input
var_names = {'${D}_{\mathrm{c}}$','$\langle\dot{\gamma}\rangle_t$'};

%% Read input file
myVars = {'INPUT','exp_design'};
load('timeCal2_00input_7000sims_19May22.mat',myVars{:})
Ns = size(exp_design,1);
M = size(exp_design,2);

%% Output
folderPath = 'TimeCalibration2_7000';
% Collect output from the solution files
fprintf('Reading the output files ... \n')
[OUTPUT,crushed_sim_idx] = timeCal_read_output(Ns,folderPath);

% delete crushed simulations and relative exp_design
OUTPUT(crushed_sim_idx) = [];
exp_design(crushed_sim_idx,:) = [];

% Transform the structure into sth readable
fprintf('Transforming output ... \n')
[outToPCE] = timeCal_getOutputReadable(OUTPUT);

% safe
fprintf('Saving the workspace ... \n')
destinationSave = 'M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_v2';
try
    cd(destinationSave)
catch
    mkdir(destinationSave)
    cd(destinationSave)
end
save('RawOutput_19May22_7000sim.mat','-v7.3')
cd(root_destination)
% %% Statistics on the output
% fprintf('Get the statistics ... \n')
% cd(root_destination)
% try
%     dest_SAplot = sprintf('Plot\\%s\\Stats',folderPath);
%     cd(dest_SAplot)
% catch
%     mkdir(dest_SAplot)
%     cd(dest_SAplot)
% end
% [timeCal_stats,outToPCE,outToPCE.diff_signal] = timeCal_statisticsOutput(outToPCE,exp_design,fittedData);
% cd(root_destination)





