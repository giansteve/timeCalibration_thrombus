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
cd M:\IFM\User\melito\Server\Projects\TimeCalibration
root_destination = pwd;
addpath(root_destination)

% initialize UQlab
% uqlab

%% Load post PCE file
load('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_validation1000\TimeCal_postSurrogate_AliModel_validation1000.mat')
M = size(INPUT.nonConst,2);
%% Surrogate accuracy display
% generate surrogate evaluations
exp_design_pce_eval = uq_getSample(INPUT,12000);
Y_pce_HS = uq_evalModel(PCE_HS,exp_design_pce_eval);
Y_pce_LS = uq_evalModel(PCE_LS,exp_design_pce_eval);
% plot
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel_validation1000\\Surrogate');
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
    dest_SAplot = sprintf('Plot_AliModel_validation1000\\SA');
    cd(dest_SAplot)
catch
    mkdir(dest_SAplot)
    cd(dest_SAplot)
end
SA_plot_time(SA_HS,linspace(0,1,size(SA_HS.main,1)),'SA_HS','$t^*$')
subplot(211); legend(var_names{INPUT.nonConst},'Interpreter','latex','Location','bestoutside')
GM_printBMP(400,400,'SA_HS')
GM_printEPS(400,400,'SA_HS')

SA_plot_time(SA_LS,linspace(0,1,size(SA_LS.main,1)),'SA_LS','$t^*$')
subplot(211); legend(var_names{INPUT.nonConst},'Interpreter','latex','Location','bestoutside')
GM_printBMP(400,400,'SA_LS')
GM_printEPS(400,400,'SA_LS')
cd(root_destination)

%% prepare for Bayesian Inverse problem
inversion_type = 'AIES';
rng(100)
clear('bayesOpts','myPriorDist','ForwardModels','myData','discrepancyOpts','myBayesian_bothModels')
clc

% 1. create PRIOR: the input PDF of my model parameters, considering the SA
% results. We want to find the POSTERIOR of only the sensitivite parameters
% bounds_prior = [min(exp_design); max(exp_design)];
% sensitivite_param = [1 6];
% nonSensitive_param = [2 3 4 5 7];
% for m = 1:M
%     if logical(sum(eq(m,nonSensitive_param)))
%         prior.Marginals(m).Name = var_names{m};
%         prior.Marginals(m).Type = 'Constant';
%         prior.Marginals(m).Parameters = mean(bounds_prior(:,m));
%     else
%         prior.Marginals(m).Name = var_names{m};
%         prior.Marginals(m).Type = 'Uniform';
%         prior.Marginals(m).Parameters = bounds_prior(:,m);
%     end
% end
% myPriorDist = uq_createInput(prior);

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
bayesOpts.Prior = INPUT;

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

if strcmpi(inversion_type,'MH')
    % 5. chose the solver
    bayesOpts.solver.Type = 'MCMC';
    bayesOpts.solver.MCMC.Sampler = 'MH'; % metropolis-hasting
    bayesOpts.solver.MCMC.Steps = 7000; % scalar to impose number of iterations
    bayesOpts.solver.MCMC.NChains = 300; % number of chains: starting point in the input domain per dimension
    % live visualization, enable only for mistuning check
    bayesOpts.solver.MCMC.Visualize.Parameters = [1;2];
    bayesOpts.solver.MCMC.Visualize.Interval = 250; % every xx steps
    % RUN IT FORREST
    myBayesian_bothModels = uq_createAnalysis(bayesOpts);
    
    % results are stored into myBayesian.Results
    % generate good posterior sample with
    uq_postProcessInversion(myBayesian_bothModels,...
        'burnIn', 0.60,... % specify the fraction of samples discarded as burn-in
        'pointEstimate',{'Mean','MAP'},... % compute 'Mean': empirical mean from sample; 'MAP': maximum posterior probability from sample
        'gelmanRubin',true... % multivariate potential scale reduction factor is computed [convergence at 1]
        )
    
    % Post-processing
    myBayesian_bothModels.Internal.FullPrior.Marginals(1).Name = 'Dc';
    %     myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = 'ct';
    myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = 'gamma';
    myBayesian_bothModels.Internal.FullPrior.Marginals(3).Name = 'e_HS';
    myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = 'e_LS'; % reChange name back to original
    uq_print(myBayesian_bothModels)
    myBayesian_bothModels.Internal.FullPrior.Marginals(1).Name = var_names{1};
    %     myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = var_names{3};
    myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = var_names{6};
    myBayesian_bothModels.Internal.FullPrior.Marginals(3).Name = '$\epsilon_{HS}$';
    myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = '$\epsilon_{LS}$'; % reChange name back to original
    
    % uq_display(myBayesian_bothModels)
    uq_display(myBayesian_bothModels,...
        'scatterplot','none',... % plot an M dimensional scatterpplot of the sample
        'trace','all',... % trace plot of MCMC chains
        'meanConvergence','all',... % convergence plot of the empirical mean
        'acceptance',true... % acceptance ratio for all chains
        )
    
    % save
    cd('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_validation1000\MH')
    save('_testMH_5k_TimeCal_postBayesian_AliModel00_gaussianDiscrepancy.mat','-v7.3')
    cd(root_destination)
    
elseif strcmpi(inversion_type,'aies') % AIES algorithm
    % 5. chose the solver
    bayesOpts.solver.Type = 'MCMC';
    bayesOpts.solver.MCMC.Sampler = 'AIES'; % metropolis-hasting
    bayesOpts.solver.MCMC.Steps = 5000; % default: 300
    bayesOpts.solver.MCMC.NChains = 200; % default: 100
    bayesOpts.solver.MCMC.a = 2; % scalar for the AIES solver
    % live visualization, enable only for mistuning check
    bayesOpts.solver.MCMC.Visualize.Parameters = [1;2];
    bayesOpts.solver.MCMC.Visualize.Interval = 500; % every xx steps
    % RUN IT FORREST
    myBayesian_bothModels = uq_createAnalysis(bayesOpts);
    
    % results are stored into myBayesian.Results
    % generate good posterior sample with
    uq_postProcessInversion(myBayesian_bothModels,...
        'burnIn', 0.60,... % specify the fraction of samples discarded as burn-in
        'pointEstimate',{'Mean','MAP'},... % compute 'Mean': empirical mean from sample; 'MAP': maximum posterior probability from sample
        'gelmanRubin',true... % multivariate potential scale reduction factor is computed [convergence at 1]
        )
    
    % Post-processing
    myBayesian_bothModels.Internal.FullPrior.Marginals(1).Name = 'Dc';
    %     myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = 'ct';
    myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = 'gamma';
    myBayesian_bothModels.Internal.FullPrior.Marginals(3).Name = 'e_HS';
    myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = 'e_LS'; % reChange name back to original
    uq_print(myBayesian_bothModels)
    myBayesian_bothModels.Internal.FullPrior.Marginals(1).Name = var_names{1};
    %     myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = var_names{3};
    myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = var_names{6};
    myBayesian_bothModels.Internal.FullPrior.Marginals(3).Name = '$\epsilon_{HS}$';
    myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = '$\epsilon_{LS}$'; % reChange name back to original
    
    % uq_display(myBayesian_bothModels)
    uq_display(myBayesian_bothModels,...
        'scatterplot','all',... % plot an M dimensional scatterpplot of the sample
        'trace','all',... % trace plot of MCMC chains
        'meanConvergence','all',... % convergence plot of the empirical mean
        'acceptance',true... % acceptance ratio for all chains
        )
    
    % save
    cd('M:\IFM\User\melito\Server\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration_validation1000\AIES')
    save('_AIES_validation1000_default_TimeCal_postBayesian_AliModel00_gaussianDiscrepancy.mat','-v7.3')
    cd(root_destination)
end




%% Plot Bayesian results
% Copying data int local variables to be used as input for
% mod_UQLab_plotSeriesPred_2021 in the following
myData_bayes = myBayesian_bothModels.Data(1);
mySamples.PostPred = myBayesian_bothModels.Results.PostProc.PostPredSample(1).PostPred;
mySamples.Post = myBayesian_bothModels.Results.PostProc.PostPredSample(1).Post;
myPointEstimate = myBayesian_bothModels.Results.PostProc.PointEstimate.ForwardRun{1,2};

PostSample3D = myBayesian_bothModels.Results.PostProc.PostSample;
PostSample2D = reshape(permute(PostSample3D, [2 1 3]), size(PostSample3D, 2), []).';

% modify limits
% PostSample2D(PostSample2D(:,1)>2e-7,:) = [];

colorRange = [.6 .6 .6;
    0.0 0.0 0.0];
n_limit = 100000;
discrepancyAsVariable = true; % if discrepancy is also inverted in Bayesian
numOutput = 2; % number of outputs

figure
if ~discrepancyAsVariable
    Y_allDim = PostSample2D(randi(size(PostSample2D,1),1,n_limit),1:end-numOutput);
    subplot_counter = 1;
    for ii = 1:size(Y_allDim,2)
        for jj = 1:size(Y_allDim,2)
            if ii == jj % HISTOGRAM PLOT
                Y = Y_allDim(:,ii);
                subplot(size(Y_allDim,2),size(Y_allDim,2),subplot_counter)
                % The width of a histogram element is computed by the Scott's rule
                w = 3.49*std(Y)*numel(Y)^(-1/3);  % Width of a histogram element
                nBins = max(ceil(range(Y)/w),1);     % Number of histograms
                [hY,hX] = hist(Y,nBins);
                % compute color for each bar
                colorFrac = hY./max(hY);
                colorVec = colorFrac'*colorRange(2,:) + (1-colorFrac')*colorRange(1,:);
                normfac = 1/(sum(hY*mean(diff(hX))));
                hY = hY*normfac;
                bar(hX,hY,'BarWidth',1,'CData',colorVec,'FaceColor','flat','EdgeAlpha',0)
                hold on
                %                 xline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1,2}(1),'r');
                xlabel(var_names{sensitivite_param(ii)})
                subplot_counter = subplot_counter + 1;
                
            elseif jj < ii % SCATTER PLOT
                Y = Y_allDim(randi(size(PostSample2D,1),1,n_limit),[jj ii]);
                subplot(size(Y_allDim,2),size(Y_allDim,2),subplot_counter)
                [N,~,~,binX,binY] = histcounts2(Y(:,1), Y(:,2));
                [NX, NY] = size(N);
                % loop over bins
                CurrColorVec = nan(size(Y,1),3);
                for kk = 1:length(N(:))
                    [xk,yk] = ind2sub(size(N),kk);
                    currFrac = N(xk,yk)/max(N(:));
                    currColor = (currFrac)*colorRange(2,:) + (1-currFrac)*colorRange(1,:);
                    currPointIds = xk == binX & yk == binY;
                    CurrColorVec(currPointIds,:) = repmat(currColor,sum(currPointIds),1);
                end
                scatter(Y(:,1), Y(:,2),2,CurrColorVec);
                % scatter(Y(:,1), Y(:,2),10,'Marker',".",'MarkerEdgeColor',[0 0 0],'MarkerEdgeAlpha',0.05);
                xlabel(var_names{sensitivite_param(jj)})
                ylabel(var_names{sensitivite_param(ii)})
                subplot_counter = subplot_counter + 1;
                
            else % leave blank subplot
                subplot_counter = subplot_counter + 1;
            end
        end
    end
else % show discrepancy inversion posterior
    Y_allDim = PostSample2D(randi(size(PostSample2D,1),1,n_limit),:);
    subplot_counter = 1;
    for ii = 1:size(Y_allDim,2)
        for jj = 1:size(Y_allDim,2)
            if ii == jj % HISTOGRAM PLOT
                Y = Y_allDim(:,ii);
                subPlotIdx_diag = diag(reshape( 1:(size(Y_allDim,2)^2),size(Y_allDim,2),size(Y_allDim,2)));
                subplot(size(Y_allDim,2),size(Y_allDim,2),subPlotIdx_diag(ii))
                % The width of a histogram element is computed by the Scott's rule
                w = 3.49*std(Y)*numel(Y)^(-1/3);  % Width of a histogram element
                nBins = max(ceil(range(Y)/w),1);     % Number of histograms
                [hY,hX] = hist(Y,nBins);
                [~,idx_max] = max(hY);
                % compute color for each bar
                colorFrac = hY./max(hY);
                colorVec = colorFrac'*colorRange(2,:) + (1-colorFrac')*colorRange(1,:);
                normfac = 1/(sum(hY*mean(diff(hX))));
                hY = hY*normfac;
                bar(hX,hY,'BarWidth',1,'CData',colorVec,'FaceColor','flat','EdgeAlpha',0)
                hold on
                maxProb_var(subplot_counter) = hX(idx_max);
                xline(maxProb_var(subplot_counter),'r','LineWidth',1);
                %                 xline(myBayesian_bothModels.Results.PostProc.PointEstimate.X{1,2}(1),'r');
                xlabel(myBayesian_bothModels.Internal.FullPrior.Marginals(ii).Name)
                subplot_counter = subplot_counter + 1;
            end
        end
    end
    for ii = 1:size(Y_allDim,2)
        for jj = 1:size(Y_allDim,2)
            if jj < ii % SCATTER PLOT
                Y = Y_allDim(:,[ii jj]);
                subPlotIdx = reshape( 1:(size(Y_allDim,2)^2),size(Y_allDim,2),size(Y_allDim,2));
                subplot(size(Y_allDim,2),size(Y_allDim,2),subPlotIdx(jj,ii))
                [N,~,~,binX,binY] = histcounts2(Y(:,1), Y(:,2));
                [NX, NY] = size(N);
                % loop over bins
                CurrColorVec = nan(size(Y,1),3);
                for kk = 1:length(N(:))
                    [xk,yk] = ind2sub(size(N),kk);
                    currFrac = N(xk,yk)/max(N(:));
                    currColor = (currFrac)*colorRange(2,:) + (1-currFrac)*colorRange(1,:);
                    currPointIds = xk == binX & yk == binY;
                    CurrColorVec(currPointIds,:) = repmat(currColor,sum(currPointIds),1);
                end
                scatter(Y(:,1), Y(:,2), 2, CurrColorVec);
                xline(maxProb_var(ii),'r','LineWidth',1);
                yline(maxProb_var(jj),'r','LineWidth',1);
                subplot_counter = subplot_counter + 1;
            else % leave blank
                subplot_counter = subplot_counter + 1;
            end
        end
    end
end

%% Inference of Posterior distribution
n_limit = 50000;
iOpts.Inference.Data = PostSample2D(randi(size(PostSample2D,1),1,n_limit),[1 2]);
iOpts.Copula.Type = 'auto';
% iOpts.Marginals(2).Type = {'Weibull','Gumbel'};
PosteriorMarginal = uq_createInput(iOpts);
posteriorSample = uq_getSample(PosteriorMarginal,10000,'lhs');
priorSample = exp_design_pce_eval(:,[1 6]);
%
figure
PosteriorData = iOpts.Inference.Data;
inferredSample = posteriorSample;
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
            xlabel(myBayesian_bothModels.Internal.FullPrior.Marginals(ii).Name)
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
            xlabel(myBayesian_bothModels.Internal.FullPrior.Marginals(ii).Name)
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
            % - Prior sample -
            subPlotIdx = reshape( 1:(size(PosteriorData,2)^2),size(PosteriorData,2),size(PosteriorData,2));
            subplot(size(PosteriorData,2),size(PosteriorData,2),subPlotIdx(jj,ii))
            scatter(Y_prior(:,1), Y_prior(:,2), 1, 'b');
            hold on
            % - Posterior sample -
            scatter(Y(:,1), Y(:,2), 1, 'r');
            % - inferred sample -
            scatter(Y_inf(:,1), Y_inf(:,2), 1, 'k');
            xlabel(myBayesian_bothModels.Internal.FullPrior.Marginals(1).Name)
            ylabel(myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name)
            subplot_counter = subplot_counter + 1;
        else % leave blank
            subplot_counter = subplot_counter + 1;
        end
    end
end


























