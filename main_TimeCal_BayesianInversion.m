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
% load('M:\IFM\User\melito\PhD\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration2_7000\TimeCal2_postSurrogate_AliModel7000.mat')
MAT = [linspace(0,1,61)' mean(phic_HS_threshold')' median(phic_HS_threshold')' mean(phic_HS_threshold')'-std(phic_HS_threshold')' mean(phic_HS_threshold')'+std(phic_HS_threshold')' mean(phic_LS_threshold')' median(phic_LS_threshold')' mean(phic_LS_threshold')'-std(phic_LS_threshold')' mean(phic_LS_threshold')'+std(phic_LS_threshold')'];
%% Surrogate accuracy display
% generate surrogate evaluations
% exp_design_pce_eval = uq_getSample(INPUT,7000,'sobol');
exp_design_pce_eval = exp_design;
Y_pce_HS = uq_evalModel(PCE_HS,exp_design_pce_eval);
Y_pce_LS = uq_evalModel(PCE_LS,exp_design_pce_eval);
% plot
cd(root_destinationC)
try
    dest_plot = sprintf('Plot_AliModel_Calibration2_timeYES\\Surrogate');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
% figure('Visible','off')
figure
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
% GM_printBMP(400,400,'ModOut_SurrOut_prob')
% GM_printEPS(400,400,'ModOut_SurrOut_prob')
% close
figure
subplot(121)
plot(phic_HS_threshold(MRI_time_index(end),:),Y_pce_HS(:,end),'k.','MarkerSize',1)
hold on
plot(linspace(0,1,100),linspace(0,1,100),'r-')
subplot(122)
plot(phic_LS_threshold(MRI_time_index(end),:),Y_pce_LS(:,end),'k.','MarkerSize',1)
hold on
plot(linspace(0,30,100),linspace(0,30,100),'r-')
cd(root_destinationC)

%% Perform sensitivity analysis
% the probability distributions of the data will not be transformed. First
% test is given in this case. In case of error, review this part.
[SA_HS.main,SA_HS.total] = SA_time(PCE_HS.PCE,M);
[SA_LS.main,SA_LS.total] = SA_time(PCE_LS.PCE,M);
cd(root_destinationC)
try
    dest_SAplot = sprintf('Plot_AliModel_Calibration2_timeYES\\SA');
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
cd(root_destinationC)

%% prepare for Bayesian Inverse problem
inversion_type = 'aies';
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
myPriorDist = INPUT;

% TRY TO COMPUTE FROM SECOND TIME STEP

% 2. give forward model, in my case the PCE surrogates
% bayesOpts_HS.ForwardModel = PCE_HS;
% bayesOpts_LS.ForwardModel = PCE_LS;
ForwardModels(1).Model = PCE_HS;
ForwardModels(1).PMap = [1 2];
ForwardModels(2).Model = PCE_LS;
ForwardModels(2).PMap = [1 2];
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

if strcmpi(inversion_type,'MH')
    % 5. chose the solver
    bayesOpts.Solver.Type = 'MCMC';
    bayesOpts.Solver.MCMC.Sampler = 'MH'; % metropolis-hasting
    bayesOpts.Solver.MCMC.Steps = 5000; % scalar to impose number of iterations
    bayesOpts.Solver.MCMC.NChains = 250; % number of chains: starting point in the input domain per dimension
    % live visualization, enable only for mistuning check
    bayesOpts.Solver.MCMC.Visualize.Parameters = [1;2];
    bayesOpts.Solver.MCMC.Visualize.Interval = 250; % every xx steps
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
    cd('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\storageFiles_noGitHub\FirstRoundCalibration\MH')
    saveName = sprintf('_MH_Steps%dNChain%d_TimeCal_AliModel00_jeffreysDiscrepancy.mat',bayesOpts.Solver.MCMC.Steps,bayesOpts.Solver.MCMC.NChains);
    save(saveName,'-v7.3')
    cd(root_destinationC)
    
elseif strcmpi(inversion_type,'aies') % AIES algorithm
    % 5. chose the solver
    bayesOpts.Solver.Type = 'MCMC';
    bayesOpts.Solver.MCMC.Sampler = 'AIES'; % metropolis-hasting
    bayesOpts.Solver.MCMC.Steps = 1000; % default: 300
    bayesOpts.Solver.MCMC.NChains = 200; % default: 100
    bayesOpts.Solver.MCMC.a = 2; % scalar for the AIES solver
    % live visualization, enable only for mistuning check
%     bayesOpts.Solver.MCMC.Visualize.Parameters = [1;2];
%     bayesOpts.Solver.MCMC.Visualize.Interval = 250; % every xx steps
    % RUN IT FORREST
    myBayesian_bothModels = uq_createAnalysis(bayesOpts);
    
    % results are stored into myBayesian.Results
    % generate good posterior sample with
    uq_postProcessInversion(myBayesian_bothModels,...
        'priorPredictive',1000,...
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
    myBayesian_bothModels.Internal.FullPrior.Marginals(2).Name = var_names{2};
    myBayesian_bothModels.Internal.FullPrior.Marginals(3).Name = '$\epsilon_{HS}$';
    myBayesian_bothModels.Internal.FullPrior.Marginals(4).Name = '$\epsilon_{LS}$'; % reChange name back to original
    
    uq_display(myBayesian_bothModels)
%     uq_display(myBayesian_bothModels,...
%         'scatterplot','all',... % plot an M dimensional scatterpplot of the sample
%         'trace','all',... % trace plot of MCMC chains
%         'meanConvergence','all',... % convergence plot of the empirical mean
%         'acceptance',true... % acceptance ratio for all chains
%         )
    
    % save
    cd('C:\Users\gm20m18\Desktop\TimeCalibrationProject\TimeCalibration\storageFiles_noGitHub\TimeCal2_Calibration_timeYES\AIES')
    saveName = sprintf('_AIES_Steps%dNChain%d_TimeCal2_AliModel00_gaussianDiscrepancy.mat',bayesOpts.Solver.MCMC.Steps,bayesOpts.Solver.MCMC.NChains);
    save(saveName,'-v7.3')
    cd(root_destinationC)
end

%% Test violin plot
% prior
hh = [];
tbl_exp_prior = [];
figure
for kk = 1:2
    subplot(1,2,kk)
    Y = myBayesian_bothModels.Results.PostProc.PriorPredSample(kk).Sample;
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
    Y = myBayesian_bothModels.Results.PostProc.PostPredSample(kk).Sample;
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
% myData_bayes = myBayesian_bothModels.Data(1);
% mySamples.PostPred = myBayesian_bothModels.Results.PostProc.PostPredSample(1).PostPred;
% mySamples.Post = myBayesian_bothModels.Results.PostProc.PostPredSample(1).Post;
% myPointEstimate = myBayesian_bothModels.Results.PostProc.PointEstimate.ForwardRun{1,2};
% 
PostSample3D = myBayesian_bothModels.Results.PostProc.PostSample;
PostSample2D = reshape(permute(PostSample3D, [2 1 3]), size(PostSample3D, 2), []).';

% modify limits
% PostSample2D(PostSample2D(:,1)>4e-7,:) = [];

colorRange = [.6 .6 .6;
    0.0 0.0 0.0];
n_limit = 5000;
discrepancyAsVariable = true; % if discrepancy is also inverted in Bayesian
numOutput = 2; % number of outputs
nonConstVec = [1,6];
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
    Y_allDim = PostSample2D(randi(size(PostSample2D,1),1,n_limit),1:end-numOutput);
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
                normfac = 1/length(Y_allDim);
                hY = hY*normfac;
                bar(hX,hY)
                hold on
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
                scatter(Y(:,1), Y(:,2), 2, [0 0 0],'MarkerEdgeAlpha',0.05);
                hold on
%                 fplot(@(x) 4.282*exp(-2.354e6*x) + 0.1245*exp(3.01e6*x))
                subplot_counter = subplot_counter + 1;
            else % leave blank
                subplot_counter = subplot_counter + 1;
            end
        end
    end
end

%% generate samples for validation
% random select element of Y_allDim
% determine how many elements is ten percent
numelements = 500;
% get the randomly-selected indices
indices = randperm(length(Y_allDim));
indices = indices(1:numelements);
% choose the subset of a you want
Y_allDim_sample = Y_allDim(indices,:);
xxx = Y_allDim_sample(:,1);
yyy = Y_allDim_sample(:,2);
figure;plot(yyy,xxx,'.');hold on;
[fitData.ff_postDcGam,fitData.gof_postDcGam] = fit(yyy,xxx,'poly2');
fplot(@(x) fitData.ff_postDcGam.p1*x^2+fitData.ff_postDcGam.p2*x+fitData.ff_postDcGam.p3,'LineWidth',1,'Color','k')
% fplot(@(x) fitData.ff_postDcGam.a*exp(fitData.ff_postDcGam.b*x) + fitData.ff_postDcGam.c*exp(fitData.ff_postDcGam.d*x),'LineWidth',1,'Color','k')
xlim([min(yyy) max(yyy)])
ylim([min(xxx) max(xxx)])


cd(root_destinationC)
save('_EDforValidationRun2_100_noInference','Y_allDim_sample')

%% save
cd('M:\IFM\User\melito\PhD\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration2_7000')
save('TimeCal2_postCalibration_preValidation.mat','-v7.3')
cd(root_destinationC)

save('_EDforValidationRun2_100_noInference','Y_allDim_sample')



















