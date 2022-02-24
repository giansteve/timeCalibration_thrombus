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
Ns = 1000; % number of simulations

% Dc
input.Marginals(1).Type = 'Gumbel';
input.Marginals(1).Parameters = [7.18e-8 2.70e-8];

% Kc
input.Marginals(2).Type = 'Constant';
input.Marginals(2).Parameters = mean([2e3 2e7]);


% Ct
input.Marginals(3).Type = 'Constant';
input.Marginals(3).Parameters = mean([1.0e+06 1.0e+07]);


% kcwall
input.Marginals(4).Type = 'Constant';
input.Marginals(4).Parameters = mean([1e2 1.00e6]);

% cbt
input.Marginals(5).Type = 'Constant';
input.Marginals(5).Parameters = mean([2e3 2.00e7]);

% gammaDot
input.Marginals(6).Type = 'Gumbel';
input.Marginals(6).Parameters = [1.022 0.2537];

% ins
input.Marginals(7).Type = 'Constant';
input.Marginals(7).Parameters = mean([0.001 5]);

INPUT = uq_createInput(input);
exp_design = uq_getSample(INPUT,Ns,'sobol');

workspace_00 = sprintf('timeCal_00input_%dsims_%s.mat',Ns,datetime('now','Format','dMMMyy'));
fprintf('Saving the workspace ... \n')
save(workspace_00,'-v7.3')

%% Create input env
model_funct_TimeCal_v2_newModel(exp_design)















