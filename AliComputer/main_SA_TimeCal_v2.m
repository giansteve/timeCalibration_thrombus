clear
close all
clc

try
    uqlab
catch
    cd /media/alireza/DATADRIVE1/Xu_Model/TimeCalibration_SA/UQLabCore_Rel1.3.0/core
    uqlab_install
    cd ../..
end

kbp = [8.000e-12 8.000e-8];
cbpt = [10 250e4];
PhiDot_c = [kbp(1)/cbpt(2) kbp(2)/cbpt(1)];
fprintf(sprintf('min(PhiDot_c): %e \t max(PhiDot_c): %e \n',min(PhiDot_c),max(PhiDot_c)))

Ns = 15; % number of simulations
%% Input

% input.Name(1) = 'Dc';
input.Marginals(1).Type = 'Uniform';
input.Marginals(1).Parameters = [1.00e-10 1.00e-6];

% input.Name(2) = 'k_c';
input.Marginals(2).Type = 'Uniform';
input.Marginals(2).Parameters = [200e2 200e4];

% input.Name(3) = 'k_BP';
input.Marginals(3).Type = 'Uniform';
input.Marginals(3).Parameters = kbp;

% input.Name(4) = 'Ct';
input.Marginals(4).Type = 'Uniform';
input.Marginals(4).Parameters = [10.00e2 10.00e4];

% input.Name(5) = 'c_BPt';
input.Marginals(5).Type = 'Uniform';
input.Marginals(5).Parameters = [20.00e2 20.00e4];

% input.Name(6) = 'T_Rt';
input.Marginals(6).Type = 'Uniform';
input.Marginals(6).Parameters = [1.00e-01 3.00e+00];

% input.Name(7) = 'k_c,wall';
input.Marginals(7).Type = 'Uniform';
input.Marginals(7).Parameters = [100 10.00e4];

% input.Name(8) = 'c_AP';
input.Marginals(8).Type = 'Uniform';
input.Marginals(8).Parameters = [1.5e14/20 4.5e14/20];

% input.Name(9) = 'c_BPbt';
input.Marginals(9).Type = 'Uniform';
input.Marginals(9).Parameters = cbpt;

% input.Name(10) = 'gammaDot_t';
input.Marginals(10).Type = 'Uniform';
input.Marginals(10).Parameters = [0.1 50];

INPUT = uq_createInput(input);
exp_design = uq_getSample(INPUT,Ns,'lhs');

% figure()
% for ii = 1:10
%     subplot(2,5,ii)
%     histogram(exp_design(:,ii),'NumBins',round(sqrt(Ns)),'Normalization','probability')
% end

workspace_00 = sprintf('timeCal_00input_%dsims_%s.mat',Ns,datetime('now','Format','dMMMyy'));
fprintf('Saving the workspace ... \n')
save(workspace_00,'-v7.3')

%% Create input env
model_funct_TimeCal_v2(exp_design)





