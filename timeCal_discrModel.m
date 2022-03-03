function [sigmaOutHS,sigmaOutLS] = timeCal_discrModel(humanData)
% create the discretization model for the ALI_newModel time calibration
% project

% SigmaOptsHS.Marginals(1).Name = 'HS_0';
% SigmaOptsHS.Marginals(1).Type = 'uniform';
% SigmaOptsHS.Marginals(1).Parameters = [humanData.H_S(1,1) 0.1];
% SigmaOptsHS.Marginals(1).Name = 'HS_1';
% SigmaOptsHS.Marginals(1).Type = 'Gaussian';
% SigmaOptsHS.Marginals(1).Parameters = [humanData.H_S(2,1) humanData.H_S(2,2)];
% SigmaOptsHS.Marginals(2).Name = 'HS_2';
% SigmaOptsHS.Marginals(2).Type = 'Gaussian';
% SigmaOptsHS.Marginals(2).Parameters = [humanData.H_S(3,1) humanData.H_S(3,2)];
% SigmaOptsHS.Marginals(3).Name = 'HS_3';
% SigmaOptsHS.Marginals(3).Type = 'Gaussian';
% SigmaOptsHS.Marginals(3).Parameters = [humanData.H_S(4,1) humanData.H_S(4,2)];
% SigmaOptsHS.Marginals(4).Name = 'HS_4';
% SigmaOptsHS.Marginals(4).Type = 'Gaussian';
% SigmaOptsHS.Marginals(4).Parameters = [humanData.H_S(5,1) humanData.H_S(5,2)];
% SigmaOptsHS.Marginals(5).Name = 'HS_5';
% SigmaOptsHS.Marginals(5).Type = 'Gaussian';
% SigmaOptsHS.Marginals(5).Parameters = [humanData.H_S(6,1) humanData.H_S(6,2)];
% SigmaOptsHS.Marginals(6).Name = 'HS_6';
% SigmaOptsHS.Marginals(6).Type = 'Gaussian';
% SigmaOptsHS.Marginals(6).Parameters = [humanData.H_S(7,1) humanData.H_S(7,2)];

SigmaOptsHS.Marginals(1).Name = '$\epsilon_{HS}$';
SigmaOptsHS.Marginals(1).Type = 'Uniform';
SigmaOptsHS.Marginals(1).Parameters = [0 0.25];

sigmaOutHS = uq_createInput(SigmaOptsHS);

% SigmaOptsLS.Marginals(1).Name = 'LS_0';
% SigmaOptsLS.Marginals(1).Type = 'uniform';
% SigmaOptsLS.Marginals(1).Parameters = [humanData.L_S(1,1) 0.1];
% SigmaOptsLS.Marginals(1).Name = 'LS_1';
% SigmaOptsLS.Marginals(1).Type = 'Gaussian';
% SigmaOptsLS.Marginals(1).Parameters = [humanData.L_S(2,1) humanData.L_S(2,2)];
% SigmaOptsLS.Marginals(2).Name = 'LS_2';
% SigmaOptsLS.Marginals(2).Type = 'Gaussian';
% SigmaOptsLS.Marginals(2).Parameters = [humanData.L_S(3,1) humanData.L_S(3,2)];
% SigmaOptsLS.Marginals(3).Name = 'LS_3';
% SigmaOptsLS.Marginals(3).Type = 'Gaussian';
% SigmaOptsLS.Marginals(3).Parameters = [humanData.L_S(4,1) humanData.L_S(4,2)];
% SigmaOptsLS.Marginals(4).Name = 'LS_4';
% SigmaOptsLS.Marginals(4).Type = 'Gaussian';
% SigmaOptsLS.Marginals(4).Parameters = [humanData.L_S(5,1) humanData.L_S(5,2)];
% SigmaOptsLS.Marginals(5).Name = 'LS_5';
% SigmaOptsLS.Marginals(5).Type = 'Gaussian';
% SigmaOptsLS.Marginals(5).Parameters = [humanData.L_S(6,1) humanData.L_S(6,2)];
% SigmaOptsLS.Marginals(6).Name = 'LS_6';
% SigmaOptsLS.Marginals(6).Type = 'Gaussian';
% SigmaOptsLS.Marginals(6).Parameters = [humanData.L_S(7,1) humanData.L_S(7,2)];

SigmaOptsLS.Marginals(1).Name = '$\epsilon_{LS}$';
SigmaOptsLS.Marginals(1).Type = 'Uniform';
SigmaOptsLS.Marginals(1).Parameters = [0 2];
sigmaOutLS = uq_createInput(SigmaOptsLS);






end 