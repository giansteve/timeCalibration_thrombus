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
cd M:\IFM\User\melito\Server\Projects\TimeCalibration
root_destination = pwd;
addpath(root_destination)

% initialize UQlab
% uqlab

%% Load 1st round
var_toLoad = 'outToPCE';
load('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_7000\AIES\_testAIES_pleaseBeTheLastOne_TimeCal_postBayesian_AliModel00_gaussianDiscrepancy.mat',var_toLoad)
MRI_time_index = round(linspace(1,61,7));
HS_1stRound(:,1) = median(outToPCE.H_S(MRI_time_index,:),2);
HS_1stRound(:,2) = std(outToPCE.H_S(MRI_time_index,:),[],2);
LS_1stRound(:,1) = median(outToPCE.L_S(MRI_time_index,:),2);
LS_1stRound(:,2) = std(outToPCE.L_S(MRI_time_index,:),[],2);
clear outToPCE

%% Load 2nd round
load('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_valid2000Copula\AIES\_validation2000Copula_done.mat')
HS_2ndRound(:,1) = mean(outToPCE.H_S(MRI_time_index,:),2);
HS_2ndRound(:,2) = std(outToPCE.H_S(MRI_time_index,:),[],2);
LS_2ndRound(:,1) = mean(outToPCE.L_S(MRI_time_index,:),2);
LS_2ndRound(:,2) = std(outToPCE.L_S(MRI_time_index,:),[],2);

%% plot fitting MRI
folderPath = 'TimeCal_MRIdata';
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel_valid2000Copula_finalization\\%s',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
xAxisShift = 0.5;
figure()
subplot(211)
hold on
% MRI data
plot(human_thr.time-xAxisShift,human_thr.H_S(:,1),'r.','MarkerSize',10)
% 1st model run
plot(human_thr.time,HS_1stRound(:,1),'kd','MarkerSize',4)
% 2nd model run
plot(human_thr.time+xAxisShift,HS_2ndRound(:,1),'kd','MarkerSize',4,'MarkerFaceColor','k')
errorbar(human_thr.time-xAxisShift,human_thr.H_S(:,1),human_thr.H_S(:,2),human_thr.H_S(:,2),'LineStyle','none','Marker','none','Color','r')
errorbar(human_thr.time,HS_1stRound(:,1),HS_1stRound(:,2),HS_1stRound(:,2),'LineStyle','none','Marker','none','Color','k')
errorbar(human_thr.time+xAxisShift,HS_2ndRound(:,1),HS_2ndRound(:,2),HS_2ndRound(:,2),'LineStyle','none','Marker','none','Color','k')
% set plot properties
ylabel('$H/S$ [-]')
% legend('MRI','$1_{st}$ model','$2_{nd}$ model','Location','se')
grid on
subplot(212)
hold on
% MRI data
plot(human_thr.time-xAxisShift,human_thr.L_S(:,1),'r.','MarkerSize',10)
% 1st model run
plot(human_thr.time,LS_1stRound(:,1),'kd','MarkerSize',4)
% 2nd model run
plot(human_thr.time+xAxisShift,LS_2ndRound(:,1),'kd','MarkerSize',4,'MarkerFaceColor','k')
errorbar(human_thr.time-xAxisShift,human_thr.L_S(:,1),human_thr.L_S(:,2),human_thr.L_S(:,2),'LineStyle','none','Marker','none','Color','r')
errorbar(human_thr.time,HS_1stRound(:,1),LS_1stRound(:,2),LS_1stRound(:,2),'LineStyle','none','Marker','none','Color','k')
errorbar(human_thr.time+xAxisShift,LS_2ndRound(:,1),LS_2ndRound(:,2),LS_2ndRound(:,2),'LineStyle','none','Marker','none','Color','k')
% posterior surrogate evaluation

% set plot properties
legend('MRI','$1_{st}$ model','$2_{nd}$ model','Location','nw')
xlabel('Time [min]')
ylabel('$L/S$ [-]')
grid on
ylim([-2 inf])
% GM_printBMP(300,300,picName)
% GM_printEPS(300,300,picName)
cd(root_destination)

%% error computations
% SSE
errRounds_HS_sse(1) = sum((human_thr.H_S(:,1) - HS_1stRound(:,1)).^2);
errRounds_HS_sse(2) = sum((human_thr.H_S(:,1) - HS_2ndRound(:,1)).^2);
errRounds_LS_sse(1) = sum((human_thr.L_S(:,1) - LS_1stRound(:,1)).^2);
errRounds_LS_sse(2) = sum((human_thr.L_S(:,1) - LS_2ndRound(:,1)).^2);
errRounds_tot_sse(1) = errRounds_HS_sse(1) + errRounds_LS_sse(1);
errRounds_tot_sse(2) = errRounds_HS_sse(2) + errRounds_LS_sse(2);

% RMS
errRounds_HS_rms(1) = sqrt(sum((human_thr.H_S(:,1) - HS_1stRound(:,1)).^2)/size(HS_1stRound,1));
errRounds_HS_rms(2) = sqrt(sum((human_thr.H_S(:,1) - HS_2ndRound(:,1)).^2)/size(HS_1stRound,1));
errRounds_LS_rms(1) = sqrt(sum((human_thr.L_S(:,1) - LS_1stRound(:,1)).^2)/size(HS_1stRound,1));
errRounds_LS_rms(2) = sqrt(sum((human_thr.L_S(:,1) - LS_2ndRound(:,1)).^2)/size(HS_1stRound,1));
errRounds_tot_rms(1) = errRounds_HS_rms(1) + errRounds_LS_rms(1);
errRounds_tot_rms(2) = errRounds_HS_rms(2) + errRounds_LS_rms(2);










