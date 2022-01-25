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
addpath('M:\IFM\User\melito\Server\Projects\safe_R1.1\safe_R1.1\PAWN')
addpath('M:\IFM\User\melito\Server\Projects\safe_R1.1\safe_R1.1\util')
addpath('M:\IFM\User\melito\Server\Projects\safe_R1.1\safe_R1.1\visualization\')

%% Script variables
excelFileName = 'WS_readDataFromMatlab_2000.xlsx';
excelFileName2 = 'WS_readDataFromMatlab.xlsx';

%% Experimental Data
% table organization
%   1   2     3      4     5     6     7     8     9
% time H/S std(H/S) L/S std(L/S) SA std(SA) Vol std(Vol)
MRI_data = readmatrix(excelFileName,'Sheet',1,'Range','B2:J9');
MRI_data(:,6:7) = MRI_data(:,6:7)./15.34e-6; % normalization of SA and vol by Step
MRI_data(:,8:9) = MRI_data(:,8:9)./15.34e-6;
%% Function fitting
% table organization
%   1   2    3   4    5
% time H/S  L/S  SA  Vol
MRI_fittedData = readmatrix(excelFileName,'Sheet',1,'Range','B12:I72');
MRI_fittedData(:,[3,5,7]) = []; % remove empty columns
MRI_fittedData(:,4) = MRI_fittedData(:,4)./15.34e-6; % normalization of SA and vol by Step
MRI_fittedData(:,5) = MRI_fittedData(:,5)./15.34e-6;
% Growth rate
MRI_growthRate = readmatrix(excelFileName,'Sheet',1,'Range','M11:P70');

%% plot fitting MRI
folderPath = 'TimeCal_Out';
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel_Calibration\\%s',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotFittingMRI(MRI_fittedData,MRI_growthRate,MRI_data,'MRI_fitting')
cd(root_destination)

%% Input
var_names = {'${D}_{\mathrm{c}}$','${k}_{\mathrm{c}}$','${k}_{\mathrm{{BP}}}$','${c}_{\mathrm{t}}$',...
    '${c}_{\mathrm{{BPt}}}$','$\overline{T}_{\mathrm{Rt}}$','${k}_{\mathrm{{c,wall}}}$',...
    '${c}_{\mathrm{{AP}}}$','${c}_{\mathrm{{BPbt}}}$','$\dot{\overline{\gamma}}_t$'};
M = size(var_names,2);
%% Read input file
% NOTE: 1500 sims and then 2000 sims
exp_design_temp2 = readmatrix(excelFileName2,'Sheet',2,'Range','B2:K3001');
[rowTrt] = find(exp_design_temp2(:,6) > 1);
exp_design_temp2(rowTrt,:) = [];
[rowGammaDot] = find(exp_design_temp2(:,10) > 20);
exp_design_temp2(rowGammaDot,:) = [];
exp_design_temp = readmatrix(excelFileName,'Sheet',2,'Range','B2:K3001');
exp_design = [exp_design_temp2;exp_design_temp];
% exp_design = exp_design_temp;
Ns = size(exp_design,1);
nonConst = 1:M;
expDesign_switch = 0;
[INPUT] = timeCal_createInput(Ns,nonConst,expDesign_switch,exp_design,0);
%% Output
Time_model_temp2 = readmatrix(excelFileName2,'Sheet',3,'Range','A3:A63');
H_S_model_temp2 = readmatrix(excelFileName2,'Sheet',3,'Range','B3:BVM63');
L_S_model_temp2 = readmatrix(excelFileName2,'Sheet',4,'Range','B3:BVM63');
Time_model_temp = readmatrix(excelFileName,'Sheet',3,'Range','A3:A63');
H_S_model_temp = readmatrix(excelFileName,'Sheet',3,'Range','B3:BVM63');
L_S_model_temp = readmatrix(excelFileName,'Sheet',4,'Range','B3:BVM63');
outToPCE.Time_model = Time_model_temp;
outToPCE.H_S_model = [H_S_model_temp2 H_S_model_temp];
outToPCE.L_S_model = [L_S_model_temp2 L_S_model_temp];
% outToPCE.H_S_model = H_S_model_temp;
% outToPCE.L_S_model = L_S_model_temp;
%% Statistics on the output
fprintf('Get the statistics ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\Stats',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
[timeCal_stats,outToPCE,outToPCE.diff_signal] = timeCal_statisticsOutput(outToPCE,exp_design,MRI_fittedData);
% scatter plot input VS likelihood
% [test,tf] = rmoutliers(outToPCE.diff_signal.likelihood_HS_LS,'mean');
figure('Visible','off')
for m = 1:M
    subplot(2,M/2,m)
%     plot(exp_design(~tf,m),(test),'k.')
    plot(exp_design(:,m),outToPCE.diff_signal.likelihood_HS_LS,'k.')
    xlabel(var_names{m})
    set(gca,'YScale','log')
    ylabel('log(L)')
end
GM_printBMP(500,300,'scatter_IN_likelihood')
GM_printEPS(500,300,'scatter_IN_likelihood')
cd(root_destination)

%% plot best model output
% this is to identify which likelihood function to use: so if to consider
% H/S in the likelihood function or not. From the results it seems that
% including H/S in the likelihood function, will drive the selection of
% many model runs which diverge in L/S.
% IN conclusion, only L/S will be considered in the LOGlikelihood function
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\Stats',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
[mval,midx] = sort(log(outToPCE.diff_signal.L_S.likelihood));
sortedED = exp_design(midx,:);
figure('Visible','off')
for cache = 1:size(outToPCE.diff_signal.L_S.likelihood,2)
    subplot(221)
    hold on
    plot(linspace(0,1,61),outToPCE.H_S_model(:,midx(cache)),'.','Color',([1 1 1] - cache/size(outToPCE.diff_signal.L_S.likelihood,2)),'MarkerSize',(8 * cache/size(outToPCE.diff_signal.L_S.likelihood,2)))
    subplot(222)
    hold on
    plot(linspace(0,1,61),outToPCE.L_S_model(:,midx(cache)),'.','Color',([1 1 1] - cache/size(outToPCE.diff_signal.L_S.likelihood,2)),'MarkerSize',(8 * cache/size(outToPCE.diff_signal.L_S.likelihood,2)))
end
subplot(221)
plot(linspace(0,1,61),MRI_fittedData(:,2),'r','LineWidth',1)
ylabel('Only L/S')
subplot(222)
plot(linspace(0,1,61),MRI_fittedData(:,3),'r','LineWidth',1)
ylim([0 10])

[mval,midx] = sort(log(outToPCE.diff_signal.likelihood_HS_LS));
sortedED = exp_design(midx,:);
for cache = 1:size(outToPCE.diff_signal.likelihood_HS_LS,2)
    subplot(223)
    hold on
    plot(linspace(0,1,61),outToPCE.H_S_model(:,midx(cache)),'.','Color',([0 0 1-cache/size(outToPCE.diff_signal.L_S.likelihood,2)]),'MarkerSize',(8 * cache/size(outToPCE.diff_signal.L_S.likelihood,2)))
    subplot(224)
    hold on
    plot(linspace(0,1,61),outToPCE.L_S_model(:,midx(cache)),'.','Color',([0 0 1-cache/size(outToPCE.diff_signal.L_S.likelihood,2)]),'MarkerSize',(8 * cache/size(outToPCE.diff_signal.L_S.likelihood,2)))
end
subplot(223)
ylabel('H/S + L/S')
plot(linspace(0,1,61),MRI_fittedData(:,2),'r','LineWidth',1)
subplot(224)
plot(linspace(0,1,61),MRI_fittedData(:,3),'r','LineWidth',1)
ylim([0 10])
GM_printBMP(400,400,'LOGlikelihood_Choice')
GM_printEPS(400,400,'LOGlikelihood_Choice')
cd(root_destination)


%% Get probability distribution of the Likelihood function output
iOpts.Inference.Data = log(outToPCE.diff_signal.L_S.likelihood)';
L_output = uq_createInput(iOpts);
uq_display(L_output)
% figure
% histogram(log(outToPCE.diff_signal.likelihood_HS_LS))

%% Metamodel
% remove outliers first, maybe also in the statistics
metamodel.Type = 'Metamodel';
metamodel.MetaType = 'PCE';
metamodel.Display = 'verbose';
% metamodel.Method = 'lars';
% metamodel.Degree = degree;
% metamodel.TruncOptions.MaxInteraction = 3;
metamodel.TruncOptions.qNorm = 0.9;
metamodel.Degree = 2:9;
metamodel.DegreeEarlyStop = false;
metamodel.Input = INPUT;
metamodel.ExpDesign.NSamples = Ns;
metamodel.ExpDesign.X = exp_design(1:Ns,:);

time_array = [1 11 21 31 41 51 61]; % indices for model instances
% Compute PCE
fprintf(' Computing PCE ... \n')

% Likelihood function
fprintf(' - LOG Likelihood L/S \n')
metamodel.ExpDesign.Y = log(outToPCE.diff_signal.L_S.likelihood)';
PCE.diff_indicat.LSlogLikelihood = uq_createModel(metamodel);

% LOG Likelihood function
% fprintf(' - LOG Likelihood \n')
% metamodel.ExpDesign.Y = log(outToPCE.diff_signal.likelihood_HS_LS)';
% PCE.diff_indicat.LOGlikelihood = uq_createModel(metamodel);

% H/S in time
fprintf(' - H/S \n')
metamodel.ExpDesign.Y = outToPCE.H_S_model(time_array,:)';
PCE.H_S.time = uq_createModel(metamodel);

% L/S in time
fprintf(' - L/S \n')
metamodel.ExpDesign.Y = outToPCE.L_S_model(time_array,:)';
PCE.L_S.time = uq_createModel(metamodel);

% % SA in time
% fprintf(' - SA \n')
% metamodel.ExpDesign.Y = outToPCE.SA_model(time_array,:)';
% PCE.SA.time = uq_createModel(metamodel);
% 
% % Vol in time
% fprintf(' - Vol \n')
% metamodel.ExpDesign.Y = outToPCE.VOL_model(time_array,:)';
% PCE.VOL.time = uq_createModel(metamodel);

% safe
fprintf('Saving the workspace ... ')
saving_time = tic;
save('timeCal_PCEdone','-v7.3')
saving_time = toc(saving_time);
fprintf(sprintf(' saving time: %.2f sec \n',saving_time))

%% Validation before calibration
% evaluate PCEs
ED_PCEinference = uq_getSample(INPUT,200000,'lhs');
PCE_infer.H_S = uq_evalModel(PCE.H_S.time,ED_PCEinference);
PCE_infer.L_S = uq_evalModel(PCE.L_S.time,ED_PCEinference);
PCE_infer.LOGlikelihood = uq_evalModel(PCE.diff_indicat.LSlogLikelihood,ED_PCEinference);
% PCE_infer.LOGlikelihood = uq_evalModel(PCE.diff_indicat.LOGlikelihood,ED_PCEinference);
% plot the results of the PCE wrt the MRI measurements and their
% correspondent error bars
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\Stats',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotPCEvsMRI(PCE_infer,outToPCE,MRI_data,'MRI_fittingPreSA')
cd(root_destination)

%% Variance-based SA on (LIKELIHOOD) and thrombus extensions in time
SA_likelihood = timeCal_SA_likelihood(PCE,M,var_names);
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\SA',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotSA_likelihood(SA_likelihood,MRI_data(:,1),M,var_names)
cd(root_destination)
close all

%% Bayes posterior updating
% What I want: update the likelihood distribution of the parameter sets
% What I need: - prior likelihood distribution of the parameter set
%              - calculated likelihood function of parameter set of new
%              observations

% [INPUT_bayes] = timeCal_createInput(Ns,nonConst,expDesign_switch,[],0);
% for cache = 2:size(time_array,2)
% % create forward models
% ModelOpts1.mFile = 'timeCal_modelComputeLikelihoodHS';
% P.MRI_fittedData = MRI_fittedData;
% P.time_instance = cache;
% P.PCE_H_S = PCE.H_S.time;
% P.PCE_L_S = PCE.L_S.time;
% P.PCE_likelihood = PCE.diff_indicat.likelihood;
% ModelOpts1.Parameters = P;
% myForwardModel(1).Model = uq_createModel(ModelOpts1);
% % myForwardModel(1).isVectorized = true;
% myForwardModel(1).PMap = 1:10; 
% %
% ModelOpts2.mFile = 'timeCal_modelComputeLikelihoodLS';
% P.MRI_fittedData = MRI_fittedData;
% P.PCE_H_S = PCE.H_S.time;
% P.PCE_L_S = PCE.L_S.time;
% ModelOpts2.Parameters = P;
% myForwardModel(2).Model = uq_createModel(ModelOpts2);
% % myForwardModel(2).isVectorized = true;
% myForwardModel(2).PMap = 1:10;
% % provide measurements
% myData(1).y = MRI_data(cache,2);
% % myData(1).y = [-10 -10];
% myData(1).Name = '$H_S$';
% myData(1).MOMap = [1;1];
% %
% myData(2).y = MRI_data(cache,4);
% myData(2).Name = '$L_S$';
% myData(2).MOMap = [2;1];
% % define discrepancy model
% DiscrepancyOpts(1).Type = 'Gaussian';
% DiscrepancyOpts(1).Parameters = MRI_data(cache,3)^2; % took from MRI excel file
% %
% DiscrepancyOpts(2).Type = 'Gaussian';
% DiscrepancyOpts(2).Parameters = MRI_data(cache,5)^2; % took from MRI excel file
% % Solver type
% Solver.Type = 'MCMC';
% Solver.MCMC.Sampler = 'AIES';
% Solver.MCMC.NChains = 40;
% Solver.MCMC.Steps = 15000;
% % Solver.MCMC.Visualize.Parameters = 1:10;
% % Solver.MCMC.Visualize.Interval = 10;
% % Implement Bayesian options
% BayesOpts.Type = 'Inversion';
% BayesOpts.Type = 'Inversion';
% BayesOpts.Data = myData;
% BayesOpts.Discrepancy = DiscrepancyOpts;
% BayesOpts.Prior = INPUT_bayes;
% BayesOpts.Data = myData;
% BayesOpts.Solver = Solver;
% BayesOpts.ForwardModel = myForwardModel;
% myBayesianAnalysis(cache) = uq_createAnalysis(BayesOpts);
% uq_print(myBayesianAnalysis(cache))
% uq_postProcessInversion(myBayesianAnalysis(cache),...
%     'posteriorPredictive',1500,...
%     'gelmanRubin',true,...          % return the convergence indicator of MCMC
%     'pointEstimate',{'Mean','MAP'}) % return the point with max posterior prob
% end
% % uq_display(myBayesianAnalysis(6))
% uq_display(myBayesianAnalysis(6),'scatterplot','all','predDist',true,'meanConvergence','all','acceptance',1)

%% try Bayes in time
% fprintf('Bayes MRI in time ... \n')
% [INPUT_bayes] = timeCal_createInput(Ns,nonConst,expDesign_switch,[],0);
% % create forward models
% ModelOpts1.mFile = 'timeCal_modelComputeLikelihoodHS';
% P.MRI_fittedData = MRI_fittedData;
% P.PCE_H_S = PCE.H_S.time;
% P.PCE_L_S = PCE.L_S.time;
% P.PCE_likelihood = PCE.diff_indicat.likelihood;
% ModelOpts1.Parameters = P;
% myForwardModel(1).Model = uq_createModel(ModelOpts1);
% myForwardModel(1).isVectorized = true;
% myForwardModel(1).PMap = 1:10; 
% %
% ModelOpts2.mFile = 'timeCal_modelComputeLikelihoodLS';
% P.MRI_fittedData = MRI_fittedData;
% P.PCE_H_S = PCE.H_S.time;
% P.PCE_L_S = PCE.L_S.time;
% ModelOpts2.Parameters = P;
% myForwardModel(2).Model = uq_createModel(ModelOpts2);
% myForwardModel(2).isVectorized = true;
% myForwardModel(2).PMap = 1:10;
% % provide measurements
% myData(1).y = MRI_data(:,2);
% myData(1).Name = '$H_S$';
% myData(1).MOMap = [1;1];
% %
% myData(2).y = MRI_data(:,4);
% myData(2).Name = '$L_S$';
% myData(2).MOMap = [2;1];
% % define discrepancy model
% DiscrepancyOpts(1).Type = 'Gaussian';
% DiscrepancyOpts(1).Parameters = mean(MRI_data(:,3).^2); % took from MRI excel file
% %
% DiscrepancyOpts(2).Type = 'Gaussian';
% DiscrepancyOpts(2).Parameters = mean(MRI_data(:,5).^2); % took from MRI excel file
% % Solver type
% Solver.Type = 'MCMC';
% Solver.MCMC.Sampler = 'AIES';
% Solver.MCMC.NChains = 30;
% Solver.MCMC.Steps = 4e4;
% % Solver.MCMC.Visualize.Parameters = 1:10;
% % Solver.MCMC.Visualize.Interval = 10;
% % Implement Bayesian options
% BayesOpts_time.Type = 'Inversion';
% BayesOpts_time.Type = 'Inversion';
% BayesOpts_time.Data = myData;
% BayesOpts_time.Discrepancy = DiscrepancyOpts;
% BayesOpts_time.Prior = INPUT_bayes;
% BayesOpts_time.Data = myData;
% BayesOpts_time.Solver = Solver;
% BayesOpts_time.ForwardModel = myForwardModel;
% myBayesianAnalysis_time = uq_createAnalysis(BayesOpts_time);
% uq_print(myBayesianAnalysis_time)
% uq_postProcessInversion(myBayesianAnalysis_time,...
%     'posteriorPredictive',1000,...
%     'gelmanRubin',true,...          % return the convergence indicator of MCMC
%     'pointEstimate',{'Mean','MAP'}) % return the point with max posterior prob
% uq_display(myBayesianAnalysis_time,'scatterplot','all','predDist',true,'meanConvergence','all','acceptance',1)
% % TEMPORARY save
% save('timeCal_Bayes_TIMEanalysis.mat','myBayesianAnalysis_time','-v7.3')
% 


%% Analyse LOG likelihood PCE evaluation
evalSample = uq_getSample(INPUT,50000,'lhs');
for m = 1:M
    if m ~= [1,3,5,6,10] % most sensitive RVs
        evalSample(:,m) = repmat(INPUT.Marginals(m).Moments(1),size(evalSample,1),1);
    end
end
% sample in good range AGAIN!!!!!!!!!!!!!
iOpts.Inference.Data = evalSample;
iOpts.Copula.Type = 'Independent';
INPUT_bayes_LLH = uq_createInput(iOpts);
% figure
% plotmatrix(evalSample)
%% Bayesian inversion on Likelihood
fprintf('Bayes LIKELIHOOD ... \n')
% [INPUT_bayes_LLH] = timeCal_createInput(Ns,nonConst,expDesign_switch,[],0);
% INPUT_bayes_LLH = INPUT;
N_mcmcSteps = 40000;
N_mcmcChain = 250;
% create forward models
ModelOpts1.mFile = 'timeCal_modelComputeLikelihood';
P.MRI_fittedData = MRI_fittedData;
P.PCE_H_S = PCE.H_S.time;
P.PCE_L_S = PCE.L_S.time;
P.PCE_likelihood = PCE.diff_indicat.LSlogLikelihood;
ModelOpts1.Parameters = P;
myForwardModel.Model = uq_createModel(ModelOpts1);
myForwardModel.PMap = 1:10; 
% provide measurements
myData(1).y = max(log(outToPCE.diff_signal.L_S.likelihood));
myData(1).Name = '$1/\sigma^2$';
% myData(1).MOMap = [1;1];
% define discrepancy model
% DiscrepancyOpts(1).Type = 'Gaussian';
% DiscrepancyOpts(1).Parameters = var(log(outToPCE.diff_signal.L_S.likelihood));
% Solver type
Solver.Type = 'MCMC';
Solver.MCMC.Sampler = 'AIES';
Solver.MCMC.NChains = N_mcmcChain;
Solver.MCMC.Steps = N_mcmcSteps;
% Solver.MCMC.Visualize.Parameters = 1:10;
% Solver.MCMC.Visualize.Interval = 10;
% Implement Bayesian options
BayesOpts_LLH.Type = 'Inversion';
BayesOpts_LLH.Data = myData;
BayesOpts_LLH.Display = 'verbose';
% BayesOpts_LLH.Discrepancy = DiscrepancyOpts;
BayesOpts_LLH.Prior = INPUT_bayes_LLH;
BayesOpts_LLH.Data = myData;
BayesOpts_LLH.Solver = Solver;
BayesOpts_LLH.ForwardModel = myForwardModel;
myBayesianAnalysis_LLH = uq_createAnalysis(BayesOpts_LLH);
uq_postProcessInversion(myBayesianAnalysis_LLH,...
    'posteriorPredictive',1000,...
    'gelmanRubin',true,...          % return the convergence indicator of MCMC
    'pointEstimate',{'Mean','MAP'}) % return the point with max posterior prob
uq_print(myBayesianAnalysis_LLH)
uq_display(myBayesianAnalysis_LLH,'scatterplot','all','trace','all','predDist',true,'meanConvergence','all','acceptance',1)
% TEMPORARY save
% save('timeCal_Bayes_TEMP_LLH_toAddToPreviousOne.mat','myBayesianAnalysis_LLH','-v7.3')


% safe
fprintf('Saving the workspace ... ')
saving_time = tic;
save('timeCal_Bayes.mat','-v7.3')
saving_time = toc(saving_time);
fprintf(sprintf(' saving time: %.2f sec \n',saving_time))

%% Get posterior distribution as if it is INPUT
% I think that to get the posterior of the input model distribution it is
% better to:
% 1. get the MAP point of the input
% 2. impose a 5% to 10% white gaussian noise to this point
% 3. produce a big data sample
% 4. get the input structure
% 5. run the PCE with this new input structure
% It should work. It is not possible to extract any information from the structure
% "myBaesianAnalysis" just because there is no information regarding the posterior.
% The only useful info is about the MAP point. So, lets use it and good night.
%%
%% sth like rand(MAP) + randn(1,1000) for the white noise
%%
%%
%% Create input from MAP
MAPpoint = myBayesianAnalysis_LLH.Results.PostProc.PointEstimate.X{2};
% generate 5% random noise from normrdn
variation1 = MAPpoint.*.01;
variation5 = MAPpoint.*.05;
variation10 = MAPpoint.*.10;
variation20 = MAPpoint.*.20;
for m = 1:size(MAPpoint,2)
MAP_1noise(:,m) = normrnd(MAPpoint(m),variation1(m),10000,1);
MAP_5noise(:,m) = normrnd(MAPpoint(m),variation5(m),10000,1);
MAP_10noise(:,m) = normrnd(MAPpoint(m),variation10(m),10000,1);
MAP_20noise(:,m) = normrnd(MAPpoint(m),variation20(m),10000,1);
end

MAP_1noise = zeros(10000,M);
cache = 1;
for m = 1:M
    if m ~= [1,3,5,6,10] % most sensitive RVs
        MAP_1noise(:,m) = repmat(INPUT.Marginals(m).Moments(1),size(evalSample,1),1);
    else
        MAP_1noise(:,m) = normrnd(MAPpoint(cache),variation1(m),10000,1);
        cache = cache + 1;
    end
end
 
% Likelihood Bayesian
iOpts.Inference.Data = MAP_1noise;
iOpts.Copula.Type = 'Independent';
PosteriorInput_LLH = uq_createInput(iOpts);
PosteriorInput_LLH = INPUT_bayes_LLH;
% uq_display(PosteriorInput_LLH)
% figure()
% histogram(myBayesianAnalysis_LLH.Results.PostProc.PostLogLikeliEval)

%% plot statistical MEAN, Median and CI of posterior RVs
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\Stats',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
for cache = 2:size(time_array,2)
    time_txt_cell{cache-1} = sprintf('MRI t=%.2f',cache/7);
end
figure('Visible','off')
for m = 1:M
    subplot(2,M/2,m)
    hold on
    mean_prior = (INPUT_bayes_LLH.Marginals(m).Parameters(1) + INPUT_bayes_LLH.Marginals(m).Parameters(2)) / 2;
    errorbar(0,mean_prior,INPUT_bayes_LLH.Marginals(m).Parameters(1)-mean_prior,INPUT_bayes_LLH.Marginals(m).Parameters(2)-mean_prior,'ro')
    % likelihood bayes
    errorbar(1,myBayesianAnalysis_LLH.Results.PostProc.Percentiles.Mean(m),...
        myBayesianAnalysis_LLH.Results.PostProc.Percentiles.Values(1,m)-myBayesianAnalysis_LLH.Results.PostProc.Percentiles.Mean(m),...
        myBayesianAnalysis_LLH.Results.PostProc.Percentiles.Values(2,m)-myBayesianAnalysis_LLH.Results.PostProc.Percentiles.Mean(m),'k.') % values from LLH
    plot(1,myBayesianAnalysis_LLH.Results.PostProc.PointEstimate.X{1,2}(m),'kx') % MAP point
    plot(1,myBayesianAnalysis_LLH.Results.PostProc.PointEstimate.X{1,1}(m),'ko') % mean point
    % MRI as vectorized output Bayes
%     errorbar(2,myBayesianAnalysis_time.Results.PostProc.Percentiles.Mean(m),...
%         myBayesianAnalysis_time.Results.PostProc.Percentiles.Values(1,m)-myBayesianAnalysis_time.Results.PostProc.Percentiles.Mean(m),...
%         myBayesianAnalysis_time.Results.PostProc.Percentiles.Values(2,m)-myBayesianAnalysis_time.Results.PostProc.Percentiles.Mean(m),'m.') % values from LLH
%     plot(2,myBayesianAnalysis_time.Results.PostProc.PointEstimate.X{1,2}(m),'mx') % MAP point
%     plot(2,myBayesianAnalysis_time.Results.PostProc.PointEstimate.X{1,1}(m),'mo') % mean point
%     xlabel(var_names{m})
    % from MRI non vectorized output Bayes
%     temp_array = [];
%     for cache = 2:size(time_array,2)
%         errorbar(cache+1,myBayesianAnalysis(cache-1).Results.PostProc.Percentiles.Mean(m),...
%             myBayesianAnalysis(cache-1).Results.PostProc.Percentiles.Values(1,m)-myBayesianAnalysis(cache-1).Results.PostProc.Percentiles.Mean(m),...
%             myBayesianAnalysis(cache-1).Results.PostProc.Percentiles.Values(2,m)-myBayesianAnalysis(cache-1).Results.PostProc.Percentiles.Mean(m),'bo') % values from MRI
%         temp_array = [temp_array; cache+1 myBayesianAnalysis(cache-1).Results.PostProc.PointEstimate.X{1,2}(m)];
% %         plot(cache,myBayesianAnalysis(cache-1).Results.PostProc.PointEstimate.X{1,2}(m),'mx','LineStyle','-') % MAP point
%     end
%     plot(temp_array(:,1),temp_array(:,2),'bx','LineStyle','-') % MAP point
%     set(gca,'XTick',0:7)
%     set(gca,'XTickLabel',['prior' 'LH' 'vecMRI' time_txt_cell])
    set(gca,'XTickLabelRotation',90)
%     if m == M
%         lgd = legend('Prior','Percen.','MAP pt','Mean pt');
%     end
end
GM_printBMP(700,500,'BayesPosterior_inputSamples')
GM_printEPS(700,500,'BayesPosterior_inputSamples')
cd(root_destination)

%% Evaluate PCE models for new posterior input
% evaluate PCEs
% posterior_ED_MRI = uq_getSample(PosteriorInput_MRI(6),3000);
posterior_ED_LLH = uq_getSample(PosteriorInput_LLH,30000,'lhs');
% posterior_ED_time = uq_getSample(PosteriorInput_time,3000);
PCE_infer.H_S = uq_evalModel(PCE.H_S.time,posterior_ED_LLH);
PCE_infer.L_S = uq_evalModel(PCE.L_S.time,posterior_ED_LLH);
% PCE_infer.SA = uq_evalModel(PCE.SA.time,posterior_ED_LLH);
% PCE_infer.VOL = uq_evalModel(PCE.VOL.time,posterior_ED_LLH);
PCE_infer.likelihood = uq_evalModel(PCE.diff_indicat.LSlogLikelihood,posterior_ED_LLH);
% plot it
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\Stats',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
figure('Visible','off')
subplot(121)
plot(linspace(0,1,61),MRI_fittedData(:,2),'r')
hold on
plot(linspace(0,1,7),mean(outToPCE.H_S_model(time_array,:),2),'k')
plot(linspace(0,1,7),mean(PCE_infer.H_S),'b')
plot(linspace(0,1,7),quantile(PCE_infer.H_S,.95),'--b')
plot(linspace(0,1,7),quantile(PCE_infer.H_S,.05),'--b')
plot(linspace(0,1,7),quantile(outToPCE.H_S_model(time_array,:),.95,2),'--k')
plot(linspace(0,1,7),quantile(outToPCE.H_S_model(time_array,:),.05,2),'--k')
lgd = legend('MRI','Model','PCE post.','Location','best');
xlabel('t/T')
ylabel('H/S')
grid on
subplot(122)
plot(linspace(0,1,61),MRI_fittedData(:,3),'r')
hold on
plot(linspace(0,1,7),mean(outToPCE.L_S_model(time_array,:),2),'k')
plot(linspace(0,1,7),mean(PCE_infer.L_S),'b')
plot(linspace(0,1,7),quantile(PCE_infer.L_S,.95),'--b')
plot(linspace(0,1,7),quantile(PCE_infer.L_S,.05),'--b')
plot(linspace(0,1,7),quantile(outToPCE.L_S_model(time_array,:),.95,2),'--k')
plot(linspace(0,1,7),quantile(outToPCE.L_S_model(time_array,:),.05,2),'--k')
grid on
xlabel('t/T')
ylabel('L/S')
GM_printBMP(300,200,'BayesPosterior_Solution')
GM_printEPS(300,200,'BayesPosterior_Solution')

% plot input posterior
figure('Visible','off')
for m = 1:M
    subplot(2,M/2,m)
    h = histcounts(posterior_ED_LLH(:,m),round(sqrt(length(posterior_ED_LLH))),'Normalization','probability');
    hold on
    plot(linspace(min(posterior_ED_LLH(:,m)),max(posterior_ED_LLH(:,m)),round(sqrt(length(posterior_ED_LLH)))),smooth(smooth(h)),'r-')
    xlabel(var_names{m})
    ylabel(sprintf('$f$(%s)',var_names{m}))
    grid on
    hold off
end
GM_printBMP(600,300,'BayesPosterior')
GM_printEPS(600,300,'BayesPosterior')

cd(root_destination)


%% Variance-based SA on (LIKELIHOOD) and thrombus extensions in time
% compute PCE of posterior first
metamodel_post.Type = 'Metamodel';
metamodel_post.MetaType = 'PCE';
metamodel_post.Display = 'verbose';
metamodel_post.Degree = 4:9;
metamodel_post.Input = PosteriorInput_LLH;
metamodel_post.ExpDesign.NSamples = size(posterior_ED_LLH,1);
metamodel_post.ExpDesign.X = posterior_ED_LLH;

% Likelihood function
metamodel_post.DegreeEarlyStop = false;
metamodel_post.ExpDesign.Y = PCE_infer.likelihood;
PCE_post.likelihood = uq_createModel(metamodel_post);

% H/S in time
metamodel_post.ExpDesign.Y = PCE_infer.H_S;
PCE_post.H_S = uq_createModel(metamodel_post);

% L/S in time
metamodel_post.ExpDesign.Y = PCE_infer.L_S;
PCE_post.L_S = uq_createModel(metamodel_post);

% Sobol inidices assessment
% H/S
[SA_post.H_S.main,SA_post.H_S.total] = SA_time_noIntegr(PCE_infer.H_S,M);
% L/S
[SA_post.L_S.main,SA_post.L_S.total] = SA_time_noIntegr(PCE_infer.L_S,M);
% Likelihood
[SA_post.likeL.main,SA_post.likeL.total,SA_post.likeL.names_sort] = SA_Sobols(PCE_infer.likelihood,M,var_names,0);

cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\SA',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotSA_likelihood(SA_post,MRI_data(:,1),M,var_names)
cd(root_destination)
close all


%% Principal Component Analysis (PCA) on posterior corr ED on LLH
% 1. standardization
posterior_ED_LLHnorm = normalize(posterior_ED_LLH);
% 2. covariance matrix
posterior_ED_LLHcov = cov(posterior_ED_LLHnorm);
posterior_ED_LLHcorr = corr(posterior_ED_LLHnorm);
% 3. eigen of covariance matrix
[V,D] = eig(posterior_ED_LLHcov); % D: eigVALUES; V: eigVECTORS
% 3.1 center the data with the new components
posterior_ED_LLHcentered = posterior_ED_LLHnorm*(V);
% 4. percentage of variance
perc_varianceLLH = diag(D) .*100 ./ sum(diag(D));
%% Principal Component Analysis (PCA) on posterior corr ED on MRI
% 1. standardization
posterior_ED_MRInorm = normalize(posterior_ED_MRI);
% 2. covariance matrix
posterior_ED_MRIcov = cov(posterior_ED_MRInorm);
posterior_ED_MRIcorr = corr(posterior_ED_MRInorm);
% 3. eigen of covariance matrix
[V,D] = eig(posterior_ED_MRIcov); % D: eigVALUES; V: eigVECTORS
% 3.1 center the data with the new components
posterior_ED_MRIcentered = posterior_ED_MRInorm*(V);
% 4. percentage of variance
perc_varianceMRI = diag(D) .*100 ./ sum(diag(D));

% plot
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\Stats',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
figure('Visible','off')
for m = 1:M
    subplot(2,M/2,m)
    plot(posterior_ED_LLHcentered(:,m),PCE_infer.likelihood,'k.')
    xlabel(var_names{m})
%     ylim([0 1])
end
GM_printBMP(400,300,'scatter_IN_PCA_likelihood')
GM_printEPS(400,300,'scatter_IN_PCA_likelihood')
figure('Visible','off')
for m = 1:M
    subplot(2,M/2,m)
    plot(posterior_ED_MRIcentered(:,m),PCE_infer.H_S(:,end),'r.')
    xlabel(var_names{m})
%     ylim([0 1])
end
GM_printBMP(400,300,'scatter_IN_PCA_HS')
GM_printEPS(400,300,'scatter_IN_PCA_HS')
figure('Visible','off')
for m = 1:M
    subplot(2,M/2,m)
    plot(posterior_ED_MRIcentered(:,m),PCE_infer.L_S(:,end),'b.')
    xlabel(var_names{m})
    ylim([0 inf])
end
GM_printBMP(400,300,'scatter_IN_PCA_LS')
GM_printEPS(400,300,'scatter_IN_PCA_LS')
cd(root_destination)





%% THE END
%% PAWN sensitivity analysis
fprintf(' Computing PAWN Sensitivity Analysis ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\SA',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
Ns_const = 150;
int_num = 15;
% Total SSE
Ereg_handle.switchMain = 0;
output_to_const = PCE.diff_indicat.totalSSE;
Y_uncond = outToPCE.diff_signal.totalSSE;
[PAWN.totSSE.Pi,PAWN.totSSE.T_m,PAWN.totSSE.T_lb,PAWN.totSSE.T_ub] = PAWNmain(Y_uncond,output_to_const,Ns_const,int_num,var_names,Ereg_handle);
pause(0.1)
grid on
ylim([0 1])
GM_printBMP(300,300,'PAWN_totSSE')
GM_printEPS(300,300,'PAWN_totSSE')

% Likelihood function
Ereg_handle.switchMain = 0;
output_to_const = PCE.diff_indicat.likelihood;
Y_uncond = outToPCE.diff_signal.likelihood_HS_LS;
[PAWN.likelihood.Pi,PAWN.likelihood.T_m,PAWN.likelihood.T_lb,PAWN.likelihood.T_ub] = PAWNmain(Y_uncond,output_to_const,Ns_const,int_num,var_names,Ereg_handle);
% pause(0.1)
grid minor
ylim([0 1])
GM_printBMP(300,300,'PAWN_likelihood')
GM_printEPS(300,300,'PAWN_likelihood')

% H/S diff indicator
Ereg_handle.switchMain = 0;
output_to_const = PCE.diff_indicat.H_S;
Y_uncond = outToPCE.diff_signal.H_S.indicator;
[PAWN.H_S_indic.Pi,PAWN.H_S_indic.T_m,PAWN.H_S_indic.T_lb,PAWN.H_S_indic.T_ub] = PAWNmain(Y_uncond,output_to_const,Ns_const,int_num,var_names,Ereg_handle);
grid on
ylim([0 1])
GM_printBMP(300,300,'PAWN_H_S_indic')
GM_printEPS(300,300,'PAWN_H_S_indic')

% L/S diff indicator
Ereg_handle.switchMain = 0;
output_to_const = PCE.diff_indicat.L_S;
Y_uncond = outToPCE.diff_signal.L_S.indicator;
[PAWN.L_S_indic.Pi,PAWN.L_S_indic.T_m,PAWN.L_S_indic.T_lb,PAWN.L_S_indic.T_ub] = PAWNmain(Y_uncond,output_to_const,Ns_const,int_num,var_names,Ereg_handle);
grid on
ylim([0 1])
GM_printBMP(300,300,'PAWN_L_S_indic')
GM_printEPS(300,300,'PAWN_L_S_indic')

% SA diff indicator
Ereg_handle.switchMain = 0;
output_to_const = PCE.diff_indicat.SA;
Y_uncond = outToPCE.diff_signal.SA.indicator;
[PAWN.SA_indic.Pi,PAWN.SA_indic.T_m,PAWN.SA_indic.T_lb,PAWN.SA_indic.T_ub] = PAWNmain(Y_uncond,output_to_const,Ns_const,int_num,var_names,Ereg_handle);
grid on
ylim([0 1])
GM_printBMP(300,300,'PAWN_SA_indic')
GM_printEPS(300,300,'PAWN_SA_indic')

% Vol diff indicator
Ereg_handle.switchMain = 0;
output_to_const = PCE.diff_indicat.Vol;
Y_uncond = outToPCE.diff_signal.Vol.indicator;
[PAWN.Vol_indic.Pi,PAWN.Vol_indic.T_m,PAWN.Vol_indic.T_lb,PAWN.Vol_indic.T_ub] = PAWNmain(Y_uncond,output_to_const,Ns_const,int_num,var_names,Ereg_handle);
grid on
ylim([0 1])
GM_printBMP(300,300,'PAWN_Vol_indic')
GM_printEPS(300,300,'PAWN_Vol_indic')
close all

cd(root_destination)
close all

%% Sensitivity analysis
SA = timeCal_SA(PCE,M,var_names);
cd(root_destination)
try
    dest_plot = sprintf('Plot\\%s\\SA',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotSA(SA,MRI_data(:,1),M,var_names)
cd(root_destination)
close all

% safe
fprintf('Saving the workspace ... ')
saving_time = tic;
save('timeCal_SAdone.mat','SA','-v7.3')
saving_time = toc(saving_time);
fprintf(sprintf(' saving time: %.2f sec \n',saving_time))

%% Evaluate PCE for Total error indicator
% input creation
% nonConst = [3,6,10];
% expDesign_switch = 1;
% [~,ED_PCEinference] = timeCal_createInput(5000,nonConst,expDesign_switch,[],0);
% Total SSE
ED_PCEinference = uq_getSample(INPUT,50000,'lhs');
PCE_infer.TotSSE = uq_evalModel(PCE.diff_indicat.totalSSE,ED_PCEinference);
[valMin_SSE,idxMin_SSE] = min(abs(PCE_infer.TotSSE));
% corresponding H/S eval
for tt = 1:length(time_array)
    H_S_PCEeval(tt) = uq_evalModel(PCE.H_S.time(tt),ED_PCEinference(idxMin_SSE,:));
    L_S_PCEeval(tt) = uq_evalModel(PCE.L_S.time(tt),ED_PCEinference(idxMin_SSE,:));
    SA_PCEeval(tt) = uq_evalModel(PCE.SA.time(tt),ED_PCEinference(idxMin_SSE,:));
    Vol_PCEeval(tt) = uq_evalModel(PCE.VOL.time(tt),ED_PCEinference(idxMin_SSE,:));
end
figure
subplot(221)
plot(MRI_data(:,2),'r:')
hold on
plot(H_S_PCEeval,'b:')
ylabel('H/S')

subplot(222)
plot(MRI_data(:,4),'r:')
hold on
plot(L_S_PCEeval,'b:')
ylabel('L/S')

subplot(223)
plot(MRI_data(:,6),'r:')
hold on
plot(SA_PCEeval,'b:')
ylabel('SA')

subplot(224)
plot(MRI_data(:,8),'r:')
hold on
plot(Vol_PCEeval,'b:')
ylabel('Vol')


% H/S indicator
ED_PCEinference = uq_getSample(INPUT,1000000,'lhs');
PCE_infer.indicH_S = uq_evalModel(PCE.diff_indicat.H_S,ED_PCEinference);
[valMin_H_S,idxMin_H_S] = min(abs(PCE_infer.indicH_S));
% corresponding H/S eval
for tt = 1:length(time_array)
    H_S_PCEeval(tt) = uq_evalModel(PCE.H_S.time(tt),ED_PCEinference(idxMin_H_S,:));
end
figure
plot(MRI_data(:,2),'r:')
hold on
plot(H_S_PCEeval,'b:')
ylabel('H/S')

% L/S indicator
ED_PCEinference = uq_getSample(INPUT,50000,'lhs');
PCE_infer.indicL_S = uq_evalModel(PCE.diff_indicat.L_S,ED_PCEinference);
[valMin_L_S,idxMin_L_S] = min(PCE_infer.indicL_S);

% SA indicator
ED_PCEinference = uq_getSample(INPUT,50000,'lhs');
PCE_infer.indicSA = uq_evalModel(PCE.diff_indicat.SA,ED_PCEinference);
[valMin_SA,idxMin_SA] = min(PCE_infer.indicSA);

% Vol indicator
ED_PCEinference = uq_getSample(INPUT,50000,'lhs');
PCE_infer.indicVol = uq_evalModel(PCE.diff_indicat.Vol,ED_PCEinference);
[valMin_Vol,idxMin_Vol] = min(PCE_infer.indicVol);


%% Reliability analysis
% H/S
Ns_RA = 5000;
fprintf('Performing Reliability Analysis ... \n')
fprintf(' - 1/4 \n')
threshold_HS = 1;
RA.H_S = timeCal_RAmain(outToPCE.diff_signal.H_S.indicator,PCE.diff_indicat.H_S,Ns_RA);

% L/S
% perform MCS RA
fprintf(' - 2/4 \n')
RA.L_S = timeCal_RAmain(outToPCE.diff_signal.L_S,PCE.diff_indicat.L_S,Ns_RA);

% SA
% perform MCS RA
fprintf(' - 3/4 \n')
RA.SA = timeCal_RAmain(outToPCE.diff_signal.SA,PCE.diff_indicat.SA,Ns_RA);

% Vol
% perform MCS RA
fprintf(' - 4/4 \n')
RA.Vol = timeCal_RAmain(outToPCE.diff_signal.Vol,PCE.diff_indicat.Vol,Ns_RA);

fprintf(' - plotting ...')
cd(root_destination)
try
    dest_RAplot = sprintf('Plot\\%s\\RA',folderPath);
    cd(dest_RAplot)
catch
    mkdir(dest_RAplot)
    cd(dest_RAplot)
end
timeCal_plotRA(RA,SA,var_names)
fprintf('done \n')
cd(root_destination)




