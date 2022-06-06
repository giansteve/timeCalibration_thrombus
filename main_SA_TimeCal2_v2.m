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
load('_EDforValidationRun2_100_noInference.mat')
%% Input
% 17/05/2022
% ----------------------------------------------
% | Parameter | par1     | par2    | type      |
% ----------------------------------------------
% | Dc        | 1e-11    | 6e-7    | Uniform   |
% | gamma     | 1e-3     | 3       | Uniform   |
% ----------------------------------------------
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

% % Dc
% input.Marginals(1).Type = 'Uniform';
% input.Marginals(1).Parameters = [1e-11 6e-7];

% % Kc
% input.Marginals(1).Type = 'Constant';
% input.Marginals(1).Parameters = mean([2e3 2e7]);
% 
% % Ct
% input.Marginals(2).Type = 'Constant';
% input.Marginals(2).Parameters = mean([1.0e+06 1.0e+07]);
% 
% % kcwall
% input.Marginals(3).Type = 'Constant';
% input.Marginals(3).Parameters = mean([1e2 1.00e6]);
% 
% % cwt
% input.Marginals(4).Type = 'Constant';
% input.Marginals(4).Parameters = mean([2e3 2.00e7]);

% % TAGammaDot
% input.Marginals(2).Type = 'Uniform';
% input.Marginals(2).Parameters = [0.001 3];

% create input of constants
% INPUT = uq_createInput(input);
exp_design = Y_allDim_sample;
Ns = length(Y_allDim_sample); % number of simulations
figure; plot(Y_allDim_sample(:,1),Y_allDim_sample(:,2),'.')
%%
workspace_00 = sprintf('timeCal2_00input_%dsims_%s.mat',Ns,datetime('now','Format','dMMMyy'));
fprintf('Saving the workspace ... \n')
save(workspace_00,'-v7.3')

%% Create input env
model_funct_TimeCal_v2_newModel(exp_design)











