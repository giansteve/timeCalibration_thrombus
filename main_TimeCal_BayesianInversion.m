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
cd C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\timeCalibration_thrombus
root_destinationC = pwd;
addpath(root_destinationC)

% initialize UQlab
% uqlab

%% Load post PCE file
load('M:\IFM\User\melito\PhD\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_7000\TimeCal_postSurrogate_AliModel7000.mat')

%% Surrogate accuracy display
% generate surrogate evaluations
% exp_design_pce_eval = uq_getSample(INPUT,7000,'sobol');
exp_design_pce_eval = exp_design;
Y_pce_HS = uq_evalModel(PCE_HS,exp_design_pce_eval);
Y_pce_LS = uq_evalModel(PCE_LS,exp_design_pce_eval);
% plot
cd(root_destinationC)
try
    dest_plot = sprintf('Plot_AliModel_Calibration\\Surrogate');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
% figure('Visible','off')
subplot(221)
GM_pdf_matrix(phic_HS_threshold(MRI_time_index,:))
xlim([0 1.5])
% ylim([0 .2])
ylabel('$p(Y)$ - model')
subplot(222)
GM_pdf_matrix(phic_LS_threshold(MRI_time_index,:))
xlim([0 35])
ylim([0 .2])
subplot(223)
GM_pdf_matrix(Y_pce_HS')
xlim([0 1.5])
% ylim([0 .2])
ylabel('$p(Y)$ - surrogate')
xlabel('$H/S$ $(-)$')
subplot(224)
GM_pdf_matrix(Y_pce_LS')
xlim([0 35])
% ylim([0 1])
xlabel('$L/S$ $(-)$')
GM_printBMP(400,400,'ModOut_SurrOut_prob')
GM_printEPS(400,400,'ModOut_SurrOut_prob')
% close
figure('Visible','off')
subplot(121)
title('$H/S$')
plot(phic_HS_threshold(MRI_time_index(end),:),Y_pce_HS(:,end),'k.','MarkerSize',1)
xlabel('Model')
ylabel('Surrogate')
hold on
plot(linspace(0,1,100),linspace(0,1,100),'r-')
subplot(122)
title('$L/S$')
plot(phic_LS_threshold(MRI_time_index(end),:),Y_pce_LS(:,end),'k.','MarkerSize',1)
hold on
xlabel('Model')
plot(linspace(0,30,100),linspace(0,30,100),'r-')
GM_printBMP(400,400,'ModOut_SurrOut')
GM_printEPS(400,400,'ModOut_SurrOut')

cd(root_destinationC)

%% Perform sensitivity analysis
% the probability distributions of the data will not be transformed. First
% test is given in this case. In case of error, review this part.
for time_inst = 1:7
    [SA_HSmain(time_inst,:),SA_HStot(time_inst,:),~] = SA_Sobols(PCE_HS.PCE(time_inst),M,var_names,0);
    [SA_LSmain(time_inst,:),SA_LStot(time_inst,:),~] = SA_Sobols(PCE_LS.PCE(time_inst),M,var_names,0);
end
cd(root_destinationC)
try
    dest_SAplot = sprintf('Plot_AliModel_Calibration\\SA');
    cd(dest_SAplot)
catch
    mkdir(dest_SAplot)
    cd(dest_SAplot)
end
figure('Visible','off')
subplot(211)
plot(SA_HSmain)
ylabel('$S_{\mathrm{{i}}}$ [-]','Interpreter','latex')
subplot(212)
plot(SA_HStot)
ylabel('$S^{\mathrm{{T}}}_{\mathrm{{i}}}$ [-]','Interpreter','latex')
xlabel('$t^*$')
% SA_plot_time(SA_HS,linspace(0,1,size(SA_HS.main,1)),'SA_HS','$t^*$')
subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
GM_printBMP(400,400,'SA_HS')
GM_printEPS(400,400,'SA_HS')

figure('Visible','off')
subplot(211)
plot(SA_LSmain)
ylabel('$S_{\mathrm{{i}}}$ [-]','Interpreter','latex')
subplot(212)
plot(SA_LStot)
ylabel('$S^{\mathrm{{T}}}_{\mathrm{{i}}}$ [-]','Interpreter','latex')
xlabel('$t^*$')
% SA_plot_time(SA_LS,linspace(0,1,size(SA_LS.main,1)),'SA_LS','$t^*$')
subplot(211); legend(var_names,'Interpreter','latex','Location','bestoutside')
GM_printBMP(400,400,'SA_LS')
GM_printEPS(400,400,'SA_LS')
cd(root_destinationC)

%% prepare for Bayesian 
NSteps = 1000;
NChain = 200;

%% prepare for Bayesian Inverse problem: HS
inversion_type = 'aies';
rng(100)
clear('bayesOpts','myPriorDist','ForwardModels','myData','discrepancyOpts','myBayesian_bothModels')
clc

% 1. create PRIOR: the input PDF of my model parameters, considering the SA
% results. We want to find the POSTERIOR of only the sensitivite parameters
% bounds_prior = [min(exp_design); max(exp_design)];
myPriorDist = INPUT;

% 4. Perform Bayesiam analysis
bayesOpts.Type = 'Inversion';
bayesOpts.Prior = myPriorDist;

% 4.1 give the discrepancy
[sigmaOptsHS,sigmaOptsLS] = timeCal_discrModel(human_thr);
discrepancyOpts(1).Type = 'Gaussian';
discrepancyOpts(1).Prior = sigmaOptsHS;
bayesOpts.Discrepancy = discrepancyOpts;

for time_inst = 1:7
    
    % 2. give forward model, in my case the PCE surrogates
    ForwardModels(1).Model = PCE_HS;
    ForwardModels(1).PMap = time_inst;
    bayesOpts.ForwardModel = ForwardModels;
    
    % 3. provide measurements
    myData(1).Name = '$H_S$';
    myData(1).y = human_thr.H_S(time_inst,1)'; % column vectors
    myData(1).MOMap = [1;...    % Model ID
        time_inst];      % Output ID
    bayesOpts.Data = myData;

    if strcmpi(inversion_type,'MH')
        % 5. chose the solver
        bayesOpts.Solver.Type = 'MCMC';
        bayesOpts.Solver.MCMC.Sampler = 'MH'; % metropolis-hasting
        bayesOpts.Solver.MCMC.Steps = 1000; % scalar to impose number of iterations
        bayesOpts.Solver.MCMC.NChains = 100; % number of chains: starting point in the input domain per dimension
        % live visualization, enable only for mistuning check
        bayesOpts.Solver.MCMC.Visualize.Parameters = [1;2];
        bayesOpts.Solver.MCMC.Visualize.Interval = 250; % every xx steps
        % RUN IT FORREST
        myBayesian_HS(time_inst) = uq_createAnalysis(bayesOpts);

        % results are stored into myBayesian.Results
        % generate good posterior sample with
        uq_postProcessInversion(myBayesian_HS(time_inst),...
            'burnIn', 0.60,... % specify the fraction of samples discarded as burn-in
            'pointEstimate',{'Mean','MAP'},... % compute 'Mean': empirical mean from sample; 'MAP': maximum posterior probability from sample
            'gelmanRubin',true... % multivariate potential scale reduction factor is computed [convergence at 1]
            )

        % Post-processing
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(1).Name = 'Dc';
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(2).Name = 'gamma';
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(3).Name = 'e_HS';
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = 'e_LS'; % reChange name back to original
        uq_print(myBayesian_HS(time_inst))
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(1).Name = var_names{1};
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = var_names{3};
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(2).Name = var_names{2};
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(3).Name = '$\epsilon_{HS}$';
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = '$\epsilon_{LS}$'; % reChange name back to original

        uq_display(myBayesian_HS(time_inst))
        %     uq_display(myBayesian_bothModels,...
        %         'scatterplot','none',... % plot an M dimensional scatterpplot of the sample
        %         'trace','all',... % trace plot of MCMC chains
        %         'meanConvergence','all',... % convergence plot of the empirical mean
        %         'acceptance',true... % acceptance ratio for all chains
        %         )

        % save
        %     cd('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\storageFiles_noGitHub\TimeCal2_Calibration\MH')
        %     saveName = sprintf('_MH_Steps%dNChain%d_TimeCal_AliModel00_jeffreysDiscrepancy.mat',bayesOpts.Solver.MCMC.Steps,bayesOpts.Solver.MCMC.NChains);
        %     save(saveName,'-v7.3')
        %     cd(root_destinationC)

    elseif strcmpi(inversion_type,'aies') % AIES algorithm
        % 5. chose the solver
        bayesOpts.Solver.Type = 'MCMC';
        bayesOpts.Solver.MCMC.Sampler = 'AIES'; % metropolis-hasting
        bayesOpts.Solver.MCMC.Steps = NSteps; % default: 300
        bayesOpts.Solver.MCMC.NChains = NChain; % default: 100
        bayesOpts.Solver.MCMC.a = 2; % scalar for the AIES solver
        % live visualization, enable only for mistuning check
%         bayesOpts.Solver.MCMC.Visualize.Parameters = [1;2];
%         bayesOpts.Solver.MCMC.Visualize.Interval = 50; % every xx steps
        % RUN IT FORREST
        myBayesian_HS(time_inst) = uq_createAnalysis(bayesOpts);

        % results are stored into myBayesian.Results
        % generate good posterior sample with
        uq_postProcessInversion(myBayesian_HS(time_inst),...
            'priorPredictive',1000,...
            'burnIn', 0.60,... % specify the fraction of samples discarded as burn-in
            'pointEstimate',{'Mean','MAP'},... % compute 'Mean': empirical mean from sample; 'MAP': maximum posterior probability from sample
            'gelmanRubin',true... % multivariate potential scale reduction factor is computed [convergence at 1]
            )

        % Post-processing
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(1).Name = 'Dc';
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(2).Name = 'gamma';
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(3).Name = 'e_HS';
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = 'e_LS'; % reChange name back to original
        uq_print(myBayesian_HS(time_inst))
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(1).Name = var_names{1};
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(2).Name = var_names{2};
        myBayesian_HS(time_inst).Internal.FullPrior.Marginals(3).Name = '$\epsilon_{HS}$';
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = '$\epsilon_{LS}$'; % reChange name back to original

%         uq_display(myBayesian_HS(time_inst))
        %     uq_display(myBayesian_bothModels,...
        %         'scatterplot','all',... % plot an M dimensional scatterpplot of the sample
        %         'trace','all',... % trace plot of MCMC chains
        %         'meanConvergence','all',... % convergence plot of the empirical mean
        %         'acceptance',true... % acceptance ratio for all chains
        %         )

        % save
        %     cd('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\storageFiles_noGitHub\TimeCal2_Calibration\AIES')
        %     saveName = sprintf('_AIES_Steps%dNChain%d_TimeCal_AliModel00_gaussianDiscrepancy.mat',bayesOpts.Solver.MCMC.Steps,bayesOpts.Solver.MCMC.NChains);
        %     save(saveName,'-v7.3')
        %     cd(root_destinationC)
    end
end

%% prepare for Bayesian Inverse problem: LS
inversion_type = 'aies';
rng(100)
clear('bayesOpts','myPriorDist','ForwardModels','myData','discrepancyOpts')
clc

% 1. create PRIOR: the input PDF of my model parameters, considering the SA
% results. We want to find the POSTERIOR of only the sensitivite parameters
% bounds_prior = [min(exp_design); max(exp_design)];
myPriorDist = INPUT;

% 4. Perform Bayesiam analysis
bayesOpts.Type = 'Inversion';
bayesOpts.Prior = myPriorDist;

% 4.1 give the discrepancy
[sigmaOptsHS,sigmaOptsLS] = timeCal_discrModel(human_thr);
discrepancyOpts(1).Type = 'Gaussian';
discrepancyOpts(1).Prior = sigmaOptsLS;
bayesOpts.Discrepancy = discrepancyOpts;

for time_inst = 1:7
    
    % 2. give forward model, in my case the PCE surrogates
    ForwardModels(1).Model = PCE_LS;
    ForwardModels(1).PMap = time_inst;
    bayesOpts.ForwardModel = ForwardModels;
    
    % 3. provide measurements
    myData(1).Name = '$L_S$';
    myData(1).y = human_thr.L_S(time_inst,1)'; % column vectors
    myData(1).MOMap = [1;...    % Model ID
        time_inst];      % Output ID
    bayesOpts.Data = myData;

    if strcmpi(inversion_type,'MH')
        % 5. chose the solver
        bayesOpts.Solver.Type = 'MCMC';
        bayesOpts.Solver.MCMC.Sampler = 'MH'; % metropolis-hasting
        bayesOpts.Solver.MCMC.Steps = 1000; % scalar to impose number of iterations
        bayesOpts.Solver.MCMC.NChains = 100; % number of chains: starting point in the input domain per dimension
        % live visualization, enable only for mistuning check
        bayesOpts.Solver.MCMC.Visualize.Parameters = [1;2];
        bayesOpts.Solver.MCMC.Visualize.Interval = 250; % every xx steps
        % RUN IT FORREST
        myBayesian_LS(time_inst) = uq_createAnalysis(bayesOpts);

        % results are stored into myBayesian.Results
        % generate good posterior sample with
        uq_postProcessInversion(myBayesian_LS(time_inst),...
            'burnIn', 0.60,... % specify the fraction of samples discarded as burn-in
            'pointEstimate',{'Mean','MAP'},... % compute 'Mean': empirical mean from sample; 'MAP': maximum posterior probability from sample
            'gelmanRubin',true... % multivariate potential scale reduction factor is computed [convergence at 1]
            )

        % Post-processing
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(1).Name = 'Dc';
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(2).Name = 'gamma';
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(3).Name = 'e_HS';
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = 'e_LS'; % reChange name back to original
        uq_print(myBayesian_LS(time_inst))
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(1).Name = var_names{1};
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = var_names{3};
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(2).Name = var_names{2};
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(3).Name = '$\epsilon_{HS}$';
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = '$\epsilon_{LS}$'; % reChange name back to original

        uq_display(myBayesian_LS(time_inst))
        %     uq_display(myBayesian_bothModels,...
        %         'scatterplot','none',... % plot an M dimensional scatterpplot of the sample
        %         'trace','all',... % trace plot of MCMC chains
        %         'meanConvergence','all',... % convergence plot of the empirical mean
        %         'acceptance',true... % acceptance ratio for all chains
        %         )

        % save
        %     cd('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\storageFiles_noGitHub\TimeCal2_Calibration\MH')
        %     saveName = sprintf('_MH_Steps%dNChain%d_TimeCal_AliModel00_jeffreysDiscrepancy.mat',bayesOpts.Solver.MCMC.Steps,bayesOpts.Solver.MCMC.NChains);
        %     save(saveName,'-v7.3')
        %     cd(root_destinationC)

    elseif strcmpi(inversion_type,'aies') % AIES algorithm
        % 5. chose the solver
        bayesOpts.Solver.Type = 'MCMC';
        bayesOpts.Solver.MCMC.Sampler = 'AIES'; % metropolis-hasting
        bayesOpts.Solver.MCMC.Steps = NSteps; % default: 300
        bayesOpts.Solver.MCMC.NChains = NChain; % default: 100
        bayesOpts.Solver.MCMC.a = 2; % scalar for the AIES solver
        % live visualization, enable only for mistuning check
%         bayesOpts.Solver.MCMC.Visualize.Parameters = [1;2];
%         bayesOpts.Solver.MCMC.Visualize.Interval = 50; % every xx steps
        % RUN IT FORREST
        myBayesian_LS(time_inst) = uq_createAnalysis(bayesOpts);

        % results are stored into myBayesian.Results
        % generate good posterior sample with
        uq_postProcessInversion(myBayesian_LS(time_inst),...
            'priorPredictive',1000,...
            'burnIn', 0.60,... % specify the fraction of samples discarded as burn-in
            'pointEstimate',{'Mean','MAP'},... % compute 'Mean': empirical mean from sample; 'MAP': maximum posterior probability from sample
            'gelmanRubin',true... % multivariate potential scale reduction factor is computed [convergence at 1]
            )

        % Post-processing
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(1).Name = 'Dc';
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(2).Name = 'gamma';
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(3).Name = 'e_HS';
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = 'e_LS'; % reChange name back to original
        uq_print(myBayesian_LS(time_inst))
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(1).Name = var_names{1};
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(2).Name = var_names{2};
        myBayesian_LS(time_inst).Internal.FullPrior.Marginals(3).Name = '$\epsilon_{HS}$';
        %     myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = '$\epsilon_{LS}$'; % reChange name back to original

%         uq_display(myBayesian_HS(time_inst))
        %     uq_display(myBayesian_bothModels,...
        %         'scatterplot','all',... % plot an M dimensional scatterpplot of the sample
        %         'trace','all',... % trace plot of MCMC chains
        %         'meanConvergence','all',... % convergence plot of the empirical mean
        %         'acceptance',true... % acceptance ratio for all chains
        %         )

        % save
        %     cd('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\storageFiles_noGitHub\TimeCal2_Calibration\AIES')
        %     saveName = sprintf('_AIES_Steps%dNChain%d_TimeCal_AliModel00_gaussianDiscrepancy.mat',bayesOpts.Solver.MCMC.Steps,bayesOpts.Solver.MCMC.NChains);
        %     save(saveName,'-v7.3')
        %     cd(root_destinationC)
    end
end

%% Test violin plot
% prior
hh = [];
tbl_exp_prior = [];
figure
for kk = 1:2
    subplot(1,2,kk)
    Y = myBayesian_HS.Results.PostProc.PriorPredSample(kk).Sample;
    for ii = 1:size(Y,2)
        % Create kernel density
        currRuns = Y(:,ii);
        [f,xi] = ksdensity(currRuns);
        % Scale and shift
        patchF = [f -fliplr(f)];
        patchF = patchF / (2*max(patchF));
        %         if exist('X','var')
        %             patchF = patchF + X(ii);
        %         else
        patchF = patchF + ii;
        %         end
        % Create the violin plot
        tbl_exp_prior = [tbl_exp_prior patchF(1:100)' flipud(patchF(101:end)') xi' fliplr(xi')];
        hh = [...
            hh;
            patch(patchF, [xi fliplr(xi)], 'b')];
    end
end
% posterior
hh = [];
tbl_exp_post = [];
hold on
for kk = 1:2
    subplot(1,2,kk)
    Y = myBayesian_HS.Results.PostProc.PostPredSample(kk).Sample;
    for ii = 1:size(Y,2)
        % Create kernel density
        currRuns = Y(:,ii);
        [f,xi] = ksdensity(currRuns);
        %         tbl_exp_post = [tbl_exp_post f' -(f') xi' fliplr(xi')];
        % Scale and shift
        patchF = [f -fliplr(f)];
        patchF = patchF / (2*max(patchF));
        %         if exist('X','var')
        %             patchF = patchF + X(ii);
        %         else
        patchF = patchF + ii;
        %         end
        % Create the violin plot
        tbl_exp_post = [tbl_exp_post patchF(1:100)' flipud(patchF(101:end)') xi' fliplr(xi')];
        hh = [...
            hh;
            patch(patchF, [xi fliplr(xi)], 'r')];
    end
end
%% Plot Bayesian results
% Copying data int local variables to be used as input for
% mod_UQLab_plotSeriesPred_2021 in the following
% myData_bayes = myBayesian_HS(1).Data(1);
% mySamples.PostPred = myBayesian_HS(1).Results.PostProc.PostPredSample(1).PostPred;
% mySamples.Post = myBayesian_HS(1).Results.PostProc.PostPredSample(1).Post;
% myPointEstimate = myBayesian_HS(1).Results.PostProc.PointEstimate.ForwardRun{1,2};

colorRange = [.6 .6 .6;
    0.0 0.0 0.0];
n_limit = 100000;
discrepancyAsVariable = true; % if discrepancy is also inverted in Bayesian
numOutput = 2; % number of outputs
nonConstVec = [1,6];

for caseNr = 1:2
    switch caseNr
        case 1 % H/S
            
            figure('Name','H/S')
            subplot_counter = 1;
            for time_inst = 1:7
                
                PostSample3D = myBayesian_HS(time_inst).Results.PostProc.PostSample;
                PostSample2D = reshape(permute(PostSample3D, [2 1 3]), size(PostSample3D, 2), []).';
                subplot(7,3,time_inst)
                
                Y_allDim = PostSample2D(randi(size(PostSample2D,1),1,n_limit),1:end-1);
                for ii = 1:size(Y_allDim,2)
                    for jj = 1:size(Y_allDim,2)
                        if ii == jj % HISTOGRAM PLOT
                            Y = Y_allDim(:,ii);
                            subPlotIdx_diag = diag(reshape( 1:(size(Y_allDim,2)^2),size(Y_allDim,2),size(Y_allDim,2)));
                            subplot(7,3,subplot_counter)
                            % The width of a histogram element is computed by the Scott's rule
                            w = 3.49*std(Y)*numel(Y)^(-1/3);  % Width of a histogram element
                            nBins = max(ceil(range(Y)/w),1);     % Number of histograms
                            [hY,hX] = hist(Y,nBins);
                            normfac = 1/(sum(hY*mean(diff(hX))));
                            hY = hY*normfac;
                            plot(hX,smooth(hY),'k')
                            hold on
                            subplot_counter = subplot_counter + 1;
                        elseif jj > ii % SCATTER PLOT
                            Y = Y_allDim(:,[ii jj]);
                            subPlotIdx = reshape( 1:(size(Y_allDim,2)^2),size(Y_allDim,2),size(Y_allDim,2));
                            subplot(7,3,subplot_counter)
                            [N,~,~,binX,binY] = histcounts2(Y(:,1), Y(:,2));
                            [NX, NY] = size(N);
                            %                             scatter(Y(:,1), Y(:,2), 2,'k','MarkerEdgeAlpha', 0.1);
                            h = histogram2(Y(:,1), Y(:,2),'DisplayStyle','tile','ShowEmptyBins','on');
                            subplot_counter = subplot_counter + 1;
                        end
                    end
                end
            end
            
        case 2 % L/S
            
            figure('Name','L/S')
            subplot_counter = 1;
            for time_inst = 1:7
                
                PostSample3D = myBayesian_LS(time_inst).Results.PostProc.PostSample;
                PostSample2D = reshape(permute(PostSample3D, [2 1 3]), size(PostSample3D, 2), []).';
                subplot(7,3,time_inst)
                
                Y_allDim = PostSample2D(randi(size(PostSample2D,1),1,n_limit),1:end-1);
                for ii = 1:size(Y_allDim,2)
                    for jj = 1:size(Y_allDim,2)
                        if ii == jj % HISTOGRAM PLOT
                            Y = Y_allDim(:,ii);
                            subPlotIdx_diag = diag(reshape( 1:(size(Y_allDim,2)^2),size(Y_allDim,2),size(Y_allDim,2)));
                            subplot(7,3,subplot_counter)
                            % The width of a histogram element is computed by the Scott's rule
                            w = 3.49*std(Y)*numel(Y)^(-1/3);  % Width of a histogram element
                            nBins = max(ceil(range(Y)/w),1);     % Number of histograms
                            [hY,hX] = hist(Y,nBins);
                            normfac = 1/(sum(hY*mean(diff(hX))));
                            hY = hY*normfac;
                            plot(hX,smooth(hY),'k')
                            subplot_counter = subplot_counter + 1;
                        elseif jj > ii % SCATTER PLOT
                            Y = Y_allDim(:,[ii jj]);
                            subPlotIdx = reshape( 1:(size(Y_allDim,2)^2),size(Y_allDim,2),size(Y_allDim,2));
                            subplot(7,3,subplot_counter)
                            [N,~,~,binX,binY] = histcounts2(Y(:,1), Y(:,2));
                            [NX, NY] = size(N);
%                             scatter(Y(:,1), Y(:,2), 2,'k','MarkerEdgeAlpha', 0.1);
                            h = histogram2(Y(:,1), Y(:,2),'DisplayStyle','tile','ShowEmptyBins','on');
                            subplot_counter = subplot_counter + 1;
                        end
                    end
                end
            end            
    end
end
%% fit posterior
[fitData.ff_postDcGam,fitData.gof_postDcGam] = fit(Y_allDim(:,1),Y_allDim(:,2),'exp2');
% random select element of Y_allDim
% determine how many elements is ten percent
numelements = 5000;
% get the randomly-selected indices
indices = randperm(length(Y_allDim));
indices = indices(1:numelements);
% choose the subset of a you want
Y_allDim_sample = Y_allDim(indices,:);

%% Inference of Posterior distribution
n_limit = 50000;
iOpts.Inference.Data = PostSample2D(randi(size(PostSample2D,1),1,n_limit),[1 2]);
iOpts.Copula.Type = 'auto';
% iOpts.Marginals(2).Type = {'Weibull','Gumbel'};
PosteriorMarginal = uq_createInput(iOpts);
posteriorSample = uq_getSample(PosteriorMarginal,5000,'lhs');
%%
% n_limit = 5000;
figure
PosteriorData = iOpts.Inference.Data;
% extractedFromPosterior = PostSample2D(randi(size(PostSample2D,1),1,n_limit),[1 2]);
inferredSample = posteriorSample;
priorSample = exp_design_pce_eval(:,[1 6]);
subplot_counter = 1;
for ii = 1:size(PosteriorData,2)
    for jj = 1:size(PosteriorData,2)
        if ii == jj % HISTOGRAM PLOT
            Y = PosteriorData(:,ii);
            Y_inf = inferredSample(:,ii);
            Y_prior = priorSample(:,ii);
            subPlotIdx_diag = diag(reshape( 1:(size(PosteriorData,2)^2),size(PosteriorData,2),size(PosteriorData,2)));
            subplot(size(PosteriorData,2),size(PosteriorData,2),subPlotIdx_diag(ii))
            % - Prior Data -
            % The width of a histogram element is computed by the Scott's rule
            w = 3.49*std(Y_prior)*numel(Y_prior)^(-1/3);  % Width of a histogram element
            nBins = max(ceil(range(Y_prior)/w),1);     % Number of histograms
            [hY,hX] = hist(Y_prior,nBins);
            [~,idx_max] = max(hY);
            normfac = 1/(sum(hY*mean(diff(hX))));
            hY = hY*normfac;
            plot(hX,hY,'b-')
            maxProb_var(subplot_counter) = hX(idx_max);
            hold on
            % - Posterior Data -
            % The width of a histogram element is computed by the Scott's rule
            w = 3.49*std(Y)*numel(Y)^(-1/3);  % Width of a histogram element
            nBins = max(ceil(range(Y)/w),1);     % Number of histograms
            [hY,hX] = hist(Y,nBins);
            [~,idx_max] = max(hY);
            normfac = 1/(sum(hY*mean(diff(hX))));
            hY = hY*normfac;
            plot(hX,hY,'r-')
            %             xline(maxProb_var(subplot_counter),'k:','LineWidth',1);
            xlabel(myBayesian_HS.Internal.FullPrior.Marginals(ii).Name)
            % - inferred sample -
            % The width of a histogram element is computed by the Scott's rule
            w = 3.49*std(Y_inf)*numel(Y_inf)^(-1/3);  % Width of a histogram element
            nBins = max(ceil(range(Y_inf)/w),1);     % Number of histograms
            [hY,hX] = hist(Y_inf,nBins);
            [~,idx_max] = max(hY);
            normfac = 1/(sum(hY*mean(diff(hX))));
            hY = hY*normfac;
            plot(hX,hY,'k')
            hold on
            maxProb_var(subplot_counter) = hX(idx_max);
            %             xline(maxProb_var(subplot_counter),'k-','LineWidth',1);
            xlabel(myBayesian_HS.Internal.FullPrior.Marginals(ii).Name)
            xlim([0 max(priorSample(:,ii))])
            legend('prior','computed posterior','inferred posterior')
            subplot_counter = subplot_counter + 1;
        end
    end
end
for ii = 1:size(PosteriorData,2)
    for jj = 1:size(PosteriorData,2)
        if jj < ii % SCATTER PLOT
            Y = PosteriorData(:,[ii jj]);
            Y_inf = inferredSample(:,[ii jj]);
            Y_prior = priorSample(:,[ii jj]);
            Y_extracted = Y_allDim_sample(:,[ii jj]);
            % - Prior sample -
            subPlotIdx = reshape( 1:(size(PosteriorData,2)^2),size(PosteriorData,2),size(PosteriorData,2));
            subplot(size(PosteriorData,2),size(PosteriorData,2),subPlotIdx(jj,ii))
            scatter(priorSample(:,1), priorSample(:,2), 1, 'b');
            hold on
            % - Posterior sample -
            scatter(PosteriorData(:,1), PosteriorData(:,2), 1, 'r');
            % - inferred sample -
            %             scatter(Y_inf(:,1), Y_inf(:,2), [], 'k');
            % - Posterior sample -
            scatter(Y_allDim_sample(:,1), Y_allDim_sample(:,2), 1, 'm');
            xlabel(myBayesian_HS.Internal.FullPrior.Marginals(1).Name)
            ylabel(myBayesian_HS.Internal.FullPrior.Marginals(2).Name)
            subplot_counter = subplot_counter + 1;
        else % leave blank
            subplot_counter = subplot_counter + 1;
        end
    end
end


%% save
cd('M:\IFM\User\melito\PhD\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_7000\AIES\04_pleaseBeTheLastOne')
save('_AIES_sim7000_1stRoundCalibration_done_newVersion.mat','-v7.3')
cd(root_destinationC)

Y_extracted = Y_allDim_sample;


save('_EDforValidationRun_5000_noInference','Y_extracted')



















