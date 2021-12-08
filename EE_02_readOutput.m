clearvars
close all
clc

set(0,'DefaultFigureWindowStyle','default')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'DefaultTextInterpreter','latex')
set(0,'defaultAxesFontSize',11)

cd M:\IFM\User\melito\Server\Projects\TimeCalibration
root_destination = pwd;
addpath(root_destination)

addpath('M:\IFM\User\melito\Server\Projects\matlab_funct\SA')
addpath('M:\IFM\User\melito\Server\Projects\matlab_funct\general')
addpath('M:\IFM\User\melito\Server\Projects\matlab_funct')

%% 3 - Model definition and solution
excelFile_input = 'EEeffectTimeCal.xlsx';
ED = readmatrix(excelFile_input,'Sheet',1,'Range','A:J');
Ns = size(ED,1);
%% Output
% cd C:\Users\gm20m18\Desktop
folderPath = 'TimeCal_EE';
% Collect output from the solution files
fprintf('Reading the output files ... \n')
[OUTPUT,crushed_sim_idx] = timeCal_read_output(Ns,folderPath); % some problem with the sims that do not run for the whole time

%% delete crushed simulations and relative exp_design
for cr_idx = 1:size(crushed_sim_idx,2)
    OUTPUT(crushed_sim_idx(cr_idx)).H_S.output = zeros(4000,1);
    OUTPUT(crushed_sim_idx(cr_idx)).L_S.output = zeros(4000,1);
end

%% Transform the structure into sth readable
fprintf('Transforming output ... \n')
thresholdThr = 0.50;
[outToPCE] = timeCal_getOutputReadable(OUTPUT,thresholdThr);

