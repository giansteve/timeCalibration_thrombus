clear
close all
clc

root_destination = pwd;
cd(root_destination)

try
    uqlab
catch
    cd /media/alireza/DATADRIVE1/Xu_Model/TimeCalibration_SA_newModel/UQLabCore_Rel1.3.0/core
    uqlab_install
    %     cd ../../
end

cd(root_destination)

kc = [8.000e-13 8.000e-11];
ct = [10 250e4];
PhiDot_c = [kc(1)/ct(2) kc(2)/ct(1)];
fprintf(sprintf('min(PhiDot_c): %e \t max(PhiDot_c): %e \n',min(PhiDot_c),max(PhiDot_c)))

%% Input
% 02/02/2022
% ---------------------------------
% | Parameter | min     | max     |
% ---------------------------------
% | Dc        | 1.0e-11 | 6.0e-07 |
% | ct        | 1.0e+06 | 1.0e+07 | 
% | gamma     | 0.001   | 3       | 
% ----------------------------------
% initial sampling of 5000
% ---------------------------------
% | Parameter | min     | max     |
% ---------------------------------
% | Dc        | 1.0e-10 | 1.0e-06 |
% | ct        | 2.0e+02 | 2.0e+06 | 
% | gamma     | 0.001   | 5       | 
% ----------------------------------
Ns = 7000; % number of simulations

% Dc
input.Marginals(1).Type = 'Uniform';
input.Marginals(1).Parameters = [1.0e-11 6.0e-07];

% Kc
input.Marginals(2).Type = 'Uniform';
input.Marginals(2).Parameters = [2e3 2e7];


% Ct
input.Marginals(3).Type = 'Uniform';
input.Marginals(3).Parameters = [1.0e+06 1.0e+07];


% kcwall
input.Marginals(4).Type = 'Uniform';
input.Marginals(4).Parameters = [1e2 1.00e6];

% cbt
input.Marginals(5).Type = 'Uniform';
input.Marginals(5).Parameters = [2e3 2.00e7];

% gammaDot
input.Marginals(6).Type = 'Uniform';
input.Marginals(6).Parameters = [0.001 3];

% ins
input.Marginals(7).Type = 'Uniform';
input.Marginals(7).Parameters = [0.001 5];

INPUT = uq_createInput(input);
exp_design = uq_getSample(INPUT,Ns,'sobol');

workspace_00 = sprintf('timeCal_00input_%dsims_%s.mat',Ns,datetime('now','Format','dMMMyy'));
fprintf('Saving the workspace ... \n')
save(workspace_00,'-v7.3')

%% Create input env
model_funct_TimeCal_v2_newModel(exp_design)


%% :::::::::::::::::::::::::: LOG RESULTS :::::::::::::::::::::::::
% Previous results 00: (02/02/2022)
%------------------- Point estimate
% --------------------------------------------------
% | Parameter | Mean    | MAP     | Parameter Type |
% --------------------------------------------------
% | Dc        | 3.8e-07 | 7.2e-08 | Model          |
% | ct        | 1.8e+06 | 1.9e+06 | Model          |
% | gamma     | 1.6     | 1.1     | Model          |
% | e_HS      | 0.031   | 0.013   | Discrepancy    |
% | e_LS      | 0.083   | 0.049   | Discrepancy    |
% --------------------------------------------------

%------------------- Correlation matrix (model parameters)
% -------------------------------------
% |       |  Dc       ct       gamma  |
% -------------------------------------
% | Dc    |  1        0.009    0.36   |
% | ct    |  0.009    1        0.4    |
% | gamma |  0.36     0.4      1      |
% -------------------------------------

%------------------- Correlation matrix (discrepancy parameters)
% -------------------------
% |      |  e_HS    e_LS  |
% -------------------------
% | e_HS |  1       0.14  |
% | e_LS |  0.14    1     |
% -------------------------














