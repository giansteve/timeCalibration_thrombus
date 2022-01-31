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
    dest_plot = sprintf('Plot_AliModel_Calibration\\%s',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotFittingMRI(fittedData,human_thr,'MRI_fitting')
cd(root_destination)

%% Define Input variables names
var_names = {'${D}_{\mathrm{c}}$',...
            '${k}_{\mathrm{c}}$',...
            '$c_\mathrm{t}$',...
            '${k}_{\mathrm{cw}}$',...
            '${c}_{\mathrm{{BPt}}}$',...
            '${\dot{\gamma}}$',...
            '${\dot{\gamma}}_\mathrm{inst}$'};

%% Read input file
myVars = {'INPUT','exp_design'};
load('timeCal_00input_5000sims_14Jan22.mat',myVars{:})

%% Load output: H/S and L/S
inp = input('Did you perform threshold optimization? [y/n] -> ','s');
if regexp(lower(inp),'y')
    load('TimeCal_OUTPUT.mat');
    clear phic_HS_threshold phic_LS_threshold OUT_time
    load('TimeCal_OUTPUT_allTime.mat')
elseif regexp(lower(inp),'n')
    load('RawOutput_25Jan22_5000sim.mat')
%     load('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\TimeCal_phic_25Jan22_5000sim.mat')
    phic_HS_threshold = outToPCE.H_S;
    phic_LS_threshold = outToPCE.L_S;
else
    error('I did not understand. \n')
end
% delete crushed sim from the input
exp_design(crushed_sim_idx,:) = [];
Ns = size(exp_design,1);
M = size(exp_design,2);

%% Statistics on the input
fprintf('Get statistics INPUT ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel_Calibration\\StatsInput');
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
    dest_plot = sprintf('Plot_AliModel_Calibration\\StatsOutput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_statisticsOutput(phic_HS_threshold,phic_LS_threshold,fittedData);
cd(root_destination)

%% Create metamodel of the thrombus model
% moving PCE computation to UQ[py]lab
MRI_time_index = round(linspace(1,61,7));
metamodel.Type = 'Metamodel';
metamodel.MetaType = 'PCE';
metamodel.Display = 'verbose';
metamodel.Method = 'lars';
metamodel.Degree = 4:10;
metamodel.TruncOptions.qNorm = [0.9 0.95];
metamodel.DegreeEarlyStop = false;
metamodel.Input = INPUT;
metamodel.ExpDesign.NSamples = Ns;
metamodel.ExpDesign.X = exp_design;
% PCE H/S
metamodel.ExpDesign.Y = phic_HS_threshold(MRI_time_index,:)';
PCE_HS = uq_createModel(metamodel);
% PCE L/S
metamodel.ExpDesign.Y = phic_LS_threshold(MRI_time_index,:)';
PCE_LS = uq_createModel(metamodel);

% save PS: remember to create function to save in specific folder with the
% date :)
save('TimeCal_postSurrogate_AliModel')

%% Surrogate accuracy display
% generate surrogate evaluations
exp_design_pce_eval = uq_getSample(INPUT,12000);
Y_pce_HS = uq_evalModel(PCE_HS,exp_design_pce_eval);
Y_pce_LS = uq_evalModel(PCE_LS,exp_design_pce_eval);
% plot
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel_Calibration\\Surrogate');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
figure('Visible','off')
subplot(221)
GM_pdf_matrix(phic_HS_threshold(MRI_time_index,:))
xlim([0 1])
% ylim([0 .2])
ylabel('$p$ $(\%)$ - model')
subplot(222)
GM_pdf_matrix(phic_LS_threshold(MRI_time_index,:))
xlim([0 inf])
% ylim([0 .2])
subplot(223)
GM_pdf_matrix(Y_pce_HS')
xlim([0 1])
% ylim([0 .2])
ylabel('$p$ $(\%)$ - surrogate')
xlabel('$H/S$ $(-)$')
subplot(224)
GM_pdf_matrix(Y_pce_LS')
xlim([0 inf])
% ylim([0 1])
xlabel('$L/S$ $(-)$')
GM_printBMP(400,400,'ModOut_SurrOut_prob')
GM_printEPS(400,400,'ModOut_SurrOut_prob')
close
cd(root_destination)

%% Perform sensitivity analysis
% the probability distributions of the data will not be transformed. First
% test is given in this case. In case of error, review this part.
[SA_HS.main,SA_HS.total] = SA_time(PCE_HS.PCE,M);
[SA_LS.main,SA_LS.total] = SA_time(PCE_LS.PCE,M);
cd(root_destination)
try
    dest_SAplot = sprintf('Plot_AliModel_Calibration\\SA');
    cd(dest_SAplot)
catch
    mkdir(dest_SAplot)
    cd(dest_SAplot)
end
SA_plot_time(SA_HS,linspace(0,1,size(SA_HS.main,1)),'SA_HS','$t^*$')
subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
GM_printBMP(400,400,'SA_HS')
GM_printEPS(400,400,'SA_HS')

SA_plot_time(SA_LS,linspace(0,1,size(SA_LS.main,1)),'SA_LS','$t^*$')
subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
GM_printBMP(400,400,'SA_LS')
GM_printEPS(400,400,'SA_LS')
cd(root_destination)

%% prepare for Bayesian Inverse problem
clear('bayesOpts','myPriorDist','ForwardModels','myData','discrepancyOpts','myBayesian_bothModels')
clc
% 1. create PRIOR: the input PDF of my model parameters, considering the SA
% results. We want to find the POSTERIOR of only the sensitivite parameters
bounds_prior = [min(exp_design); max(exp_design)];
nonSensitive_param = [2 4 5 7];
for m = 1:M
    if logical(sum(eq(m,nonSensitive_param)))
        prior.Marginals(m).Name = var_names{m};
        prior.Marginals(m).Type = 'Constant';
        prior.Marginals(m).Parameters = mean(bounds_prior(:,m));
    else
        prior.Marginals(m).Name = var_names{m};
        prior.Marginals(m).Type = 'Uniform';
        prior.Marginals(m).Parameters = bounds_prior(:,m);
    end
end
myPriorDist = uq_createInput(prior);

% TRY TO COMPUTE FROM SECOND TIME STEP

% 2. give forward model, in my case the PCE surrogates
% bayesOpts_HS.ForwardModel = PCE_HS;
% bayesOpts_LS.ForwardModel = PCE_LS;
ForwardModels(1).Model = PCE_HS;
ForwardModels(1).PMap = [1 2 3 4 5 6 7];
ForwardModels(2).Model = PCE_LS;
ForwardModels(2).PMap = [1 2 3 4 5 6 7];
bayesOpts.ForwardModel = ForwardModels;

% 3. provide measurements
% myData_HS.Name = 'H_S';
% myData_HS.y = human_thr.H_S(:,1); % column vectors
% myData_LS.Name = 'L_S';
% myData_LS.y = human_thr.L_S(:,1); % column vectors
myData(1).Name = '$H_S$';
myData(1).y = human_thr.H_S(:,1)'; % column vectors
myData(1).MOMap = [1 1 1 1 1 1 1;...    % Model ID
                   1 2 3 4 5 6 7];      % Output ID
myData(2).Name = '$L_S$';
myData(2).y = human_thr.L_S(:,1)'; % column vectors
myData(2).MOMap = [2 2 2 2 2 2 2;...    % Model ID
                   1 2 3 4 5 6 7];      % Output ID
               
% 4. Perform Bayesiam analysis
% bayesOpts_HS.Type = 'Inversion';
% bayesOpts_HS.Data = myData_HS;
% bayesOpts_LS.Type = 'Inversion';
% bayesOpts_LS.Data = myData_LS;
bayesOpts.Type = 'Inversion';
bayesOpts.Data = myData;
bayesOpts.Prior = myPriorDist;

% 4.1 give the discrepancy
% :::::::::::::::: ORIGINAL :::::::::::::::::::::::
% Discrepancy is the STD of the data. The discrepancy simulates the measurement error!
% In case of multiple outputs, the discrepancy has the shape of a vector
% whose dimension is equal to the number of dimensions. In this case: 1x7
% :::::::::::::::: 28/01/22 UPDATE ::::::::::::::::
% the discrepancy model has known mean and std at each time step. So we
% need to build it accordingly.
% :::::::::::::::: 28/01/22 UPDATE 2 ::::::::::::::
% The discrepancy is the error between the computational model and the
% data.
% % discrepancyOpts(1).Type = 'Gaussian';
% % discrepancyOpts(1).Parameters = (human_thr.H_S(:,2)').^2; % remember to square it
% % discrepancyOpts(2).Type = 'Gaussian';
% % discrepancyOpts(2).Parameters = (human_thr.L_S(:,2)').^2; % remember to square it
[sigmaOptsHS,sigmaOptsLS] = timeCal_discrModel(human_thr);
discrepancyOpts(1).Type = 'Gaussian';
discrepancyOpts(1).Prior = sigmaOptsHS;
discrepancyOpts(2).Type = 'Gaussian';
discrepancyOpts(2).Prior = sigmaOptsLS;
bayesOpts.Discrepancy = discrepancyOpts;

% chose the solver
bayesOpts.solver.Type = 'MCMC';
bayesOpts.solver.MCMC.Sampler = 'MH'; % metropolis-hasting
bayesOpts.solver.MCMC.Steps = 500000; % scalar to impose number of iterations
bayesOpts.solver.MCMC.NChains = 250; % number of chains: starting point in the input domain per dimension
% live visualization, enable only for mistuning check
bayesOpts.solver.MCMC.Visualize.Parameters = [1;2;3];
bayesOpts.solver.MCMC.Visualize.Interval = 250; % every xx steps
% RUN IT FORREST
myBayesian_bothModels = uq_createAnalysis(bayesOpts);

% results are stored into myBayesian.Results
% generate good posterior sample with
uq_postProcessInversion(myBayesian_bothModels,'burnIn', 0.8,'pointEstimate',{'Mean','MAP'})
% Post-processing
uq_print(myBayesian_bothModels)
uq_display(myBayesian_bothModels)

cd('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles')
save('_testBig_TimeCal_postBayesian_AliModel00_gaussianDiscrepancy.mat','-v7.3')
cd(root_destination)

%% Plot Bayesian results
% Copying data int local variables to be used as input for
% mod_UQLab_plotSeriesPred_2021 in the following
myData_bayes = myBayesian_bothModels.Data(1);
mySamples.PostPred = myBayesian_bothModels.Results.PostProc.PostPredSample(1).PostPred;
mySamples.Post = myBayesian_bothModels.Results.PostProc.PostPredSample(1).Post;
myPointEstimate = myBayesian_bothModels.Results.PostProc.PointEstimate.ForwardRun{1,2};

PostSample3D = myBayesian_bothModels.Results.PostProc.PostSample;
PostSample2D = reshape(permute(PostSample3D, [2 1 3]), size(PostSample3D, 2), []).';

n_limit = 1000000;
figure
subplot(221)
histogram(PostSample2D(end-n_limit:end,1),'Normalization','probability')
hold on
xline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1,2}(1),'r');
xlabel(var_names{1})
subplot(224)
histogram(PostSample2D(end-n_limit:end,2),'Normalization','probability')
hold on
xline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1,2}(2),'r');
xlabel(var_names{2})
camroll(90)
subplot(223)
scatter(PostSample2D(end-n_limit:end,1),PostSample2D(end-n_limit:end,2),1,'k')
hold on
xline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1,2}(1),'r');
yline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1,2}(2),'r');
xlabel(var_names{1})
ylabel(var_names{2})


figure
subplot(221)
histogram(PostSample2D(:,1),'Normalization','probability')
hold on
xline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1}(1),'r');
subplot(224)
histogram(PostSample2D(:,2),'Normalization','probability')
hold on
xline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1}(2),'r');
camroll(-90)
subplot(223)
scatter(PostSample2D(end-n_limit:end,1),PostSample2D(end-n_limit:end,2),1,'k')
hold on
xline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1}(1),'r');
yline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1}(2),'r');
xlabel(var_names{1})
ylabel(var_names{2})












