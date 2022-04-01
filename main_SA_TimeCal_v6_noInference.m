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

% load WS
load('_EDforValidationRun_5000_noInference.mat')
%% Input
% 20/03/2022
% ----------------------------------------------
% | Parameter | type                           |
% ----------------------------------------------
% | Dc        | extracted from posterior       |
% | gamma     | extracted from posterior       |
% ----------------------------------------------
% 02/03/2022
% ----------------------------------------------
% | Parameter | par1     | par2    | type      |
% ----------------------------------------------
% | Dc        | 7.18e-8  | 2.70e-8 | Gumbel    |
% | gamma     | 1.022    | 0.2537  | Gumbel    |
% ----------------------------------------------
% | Copula    |          |         | Pair      |
% ----------------------------------------------
% 24/02/2022
% ------------------------------------------
% | Parameter | par1    | par2    | type   |
% ------------------------------------------
% | Dc        | 7.18e-8 | 2.70e-8 | Gumbel |
% | gamma     | 1.022   | 0.2537  | Gumbel |
% ------------------------------------------
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
Ns = 5000; % number of simulations

% INPUT = constant variables;

% Kc
input.Marginals(1).Type = 'Constant';
input.Marginals(1).Parameters = mean([2e3 2e7]);

% Ct
input.Marginals(2).Type = 'Constant';
input.Marginals(2).Parameters = mean([1.0e+06 1.0e+07]);

% kcwall
input.Marginals(3).Type = 'Constant';
input.Marginals(3).Parameters = mean([1e2 1.00e6]);

% cbt
input.Marginals(4).Type = 'Constant';
input.Marginals(4).Parameters = mean([2e3 2.00e7]);

% ins
input.Marginals(5).Type = 'Constant';
input.Marginals(5).Parameters = mean([0.001 5]);

% create input of constants
INPUT = uq_createInput(input);
exp_design_const = uq_getSample(INPUT,Ns,'sobol');

% INPUT = PosteriorMarginal;
exp_design_var = Y_extracted;

% merge
exp_design = [exp_design_var(:,1) exp_design_const(:,(1:4)) exp_design_var(:,2) exp_design_const(:,5)];

workspace_00 = sprintf('timeCal_00input_%dsims_%s.mat',Ns,datetime('now','Format','dMMMyy'));
fprintf('Saving the workspace ... \n')
save(workspace_00,'-v7.3')

%% Create input env
model_funct_TimeCal_v2_newModel(exp_design)











