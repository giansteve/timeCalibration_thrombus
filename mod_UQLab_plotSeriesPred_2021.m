function mod_UQLab_plotSeriesPred_2021(Sample, Data, plotType, varargin)
% function plotSeriesPred from UQLab File
%         uq_display_uq_inversion.m (v.1.3)
% with some modifications  by O. Klein, WIAS, Berlin
%% ------------------------------------------------------------------------
%PLOTSERIESPRED creates a plot of serial (non-scalar) predictions.
% can also be used for scalar predictions 

%% Parse inputs

StartInd = 1;
EndInd = size(Data.y,2);
if nargin == 3
    % Point estimate and Index boundaries not available
    pointEstimatePredFlag = false;
    
else
    % Point estimate available
    pointEstimatePredFlag = true;
    pointEstimatePred = varargin{1};
    if nargin > 4
        StartInd =  max(1,varargin{2});
        EndInd   = min(EndInd,varargin{3});
        if StartInd > EndInd
            fprintf('mod_uq_plotSeriesPred_1: StartInd = %g > EndInd = %g \n',...
                StartInd,EndInd);
            error('This relationship is not allowed');
        end % if
    end % else
end % else


% Set custom colors
priorColor = [0.5, 0.7, 1.0];
postColor  = [0.0, 0.2, 0.6];
postModColor = [1, 170/255, 1];
postModPEplusColor = [0.8, 0.7, 0.5];
%% Create the plot
uq_figure  % Open a new figure

% Use violin plot for hist/ogram of data series
switch plotType
    
    case 'prior&Post'
        % Both prior and posterior predictive runs
        priorRuns = Sample.PriorPred(:,StartInd:EndInd);
        postRuns = Sample.PostPred(:,StartInd:EndInd);
        % Plot prior predictive runs
        priorPlot = uq_violinplot(StartInd:EndInd,priorRuns, 'FaceColor', priorColor);
        hold on
        % Plot posterior predictive runs
        postPlot = uq_violinplot(StartInd:EndInd,postRuns, 'FaceColor', postColor);
        % Update legend target
        legendObj = [priorPlot(1) postPlot(1)];
        legendName = {'prior predictive', 'posterior predictive'};
    case 'prior'
        % Only prior prior predictive runs
        priorRuns = Sample.PriorPred(:,StartInd:EndInd);
        % Plot prior predictive runs
        priorPlot = uq_violinplot(StartInd:EndInd,priorRuns, 'FaceColor', priorColor);
        % Update legend target
        legendObj = priorPlot(1);
        legendName = {'prior predictive'};
        hold on
    case 'post'
        % Only posterior predictive runs
        postRuns = Sample.PostPred(:,StartInd:EndInd);
        % Plot posterior predictive runs
        postPlot = uq_violinplot(StartInd:EndInd,postRuns, 'FaceColor', postColor);
        % Update legend target
        legendObj = postPlot(1);
        legendName = {'posterior predictive'};
        hold on
    case 'post&Model'
        % Both posterior predictive and model at samples runs
        postRuns = Sample.PostPred(:,StartInd:EndInd);
        ModelRuns = Sample.Post(:,StartInd:EndInd);
        
        % Plot posterior predictive runs combined with density for model at
        % samples
        postPlot = uq_violinplot(StartInd:EndInd,postRuns, 'FaceColor', postColor);
        hold on
        % plot model output at sample
        ModelPlot = uq_violinplot(StartInd:EndInd,ModelRuns, 'FaceColor', postModColor);
        hold on
        
        % Update legend target
        legendObj = [postPlot(1) ModelPlot(1)];
        legendName = {'posterior predictive','model at posterior'};
    case 'model'
        % Both posterior predictive and model at samples runs
        ModelRuns = Sample.Post(:,StartInd:EndInd);
        
        
        % plot model output at sample
        ModelPlot = uq_violinplot(StartInd:EndInd,ModelRuns, 'FaceColor', postModColor);
        hold on
        
        % Update legend target
        legendObj = [ModelPlot(1)];
        legendName = {'model at posterior'};
        
        
    case 'Model(PE)+noise'
        % Both posterior predictive and model at samples runs
        if pointEstimatePredFlag
            % assign model runs to current Data based on MOMap
            pointEstimatePredCurr = zeros(size(pointEstimatePred{1}(1).Out,1),size(Data.y,2));
            % loop over data points
            for ii = 1:size(pointEstimatePredCurr,2)
                currModel = Data.MOMap(1,ii);
                currOut = Data.MOMap(2,ii);
                pointEstimatePredCurr(:,ii) = pointEstimatePred{1}(currModel).Out(:,currOut);
            end
            ModelRuns = Sample.Post(:,StartInd:EndInd);
            postRuns =  Sample.PostPred(:,StartInd:EndInd);
            usedNoise = postRuns -ModelRuns;
            Model_at_PE_plus_noise = pointEstimatePredCurr(1,StartInd:EndInd) + usedNoise;
            
            % plot model output at sample
            Model_at_PE_plusPlot = uq_violinplot(StartInd:EndInd,...
                Model_at_PE_plus_noise, 'FaceColor',postModPEplusColor);
            hold on
            
            % Update legend target
            legendObj = [Model_at_PE_plusPlot(1)];
            legendName = { ['(model at ',lower(pointEstimatePred{1}(1).Type),') +noise']};
        else
            close;
            error('mod_UQLab_plotSeriesPred: Plot option "Model(PE)+noise" can only be used if point estimate is known');
        end % else to  if    if pointEstimatePredFlag
    case 'post&Model(PE)+noise'
        % Plotposterior predictive  density for (model at point estimate ) + noise
        if pointEstimatePredFlag
            % assign model runs to current Data based on MOMap
            pointEstimatePredCurr = zeros(size(pointEstimatePred{1}(1).Out,1),size(Data.y,2));
            % loop over data points
            for ii = 1:size(pointEstimatePredCurr,2)
                currModel = Data.MOMap(1,ii);
                currOut = Data.MOMap(2,ii);
                pointEstimatePredCurr(:,ii) = pointEstimatePred{1}(currModel).Out(:,currOut);
            end
            ModelRuns = Sample.Post(:,StartInd:EndInd);
            postRuns =  Sample.PostPred(:,StartInd:EndInd);
            
            Model_at_PE_plus_noise = pointEstimatePredCurr(1,StartInd:EndInd) + postRuns - ModelRuns;
            
            % plot model output at sample
            postPlot = uq_violinplot(StartInd:EndInd,postRuns, 'FaceColor', postColor);
            
            hold on
            Model_at_PE_plusPlot = uq_violinplot(StartInd:EndInd,...
                Model_at_PE_plus_noise,'FaceColor',postModPEplusColor);
            
            % Update legend target
            legendObj = [postPlot(1) Model_at_PE_plusPlot(1)];
            legendName = {'posterior predictive',...
                ['(model at ',lower(pointEstimatePred{1}(1).Type),') +noise']};
            
        else
            error('mod_UQLab_plotSeriesPred: Plot option "post&Model(PE)+noise" can only be used if point estimate is known');
        end % else to  if    if pointEstimatePredFlag
    case 'post&Model(PE)+noise:edges'
        % Plotposterior predictive  density for (model at point estimate ) + noise
        if pointEstimatePredFlag
            % assign model runs to current Data based on MOMap
            pointEstimatePredCurr = zeros(size(pointEstimatePred{1}(1).Out,1),size(Data.y,2));
            % loop over data points
            for ii = 1:size(pointEstimatePredCurr,2)
                currModel = Data.MOMap(1,ii);
                currOut = Data.MOMap(2,ii);
                pointEstimatePredCurr(:,ii) = pointEstimatePred{1}(currModel).Out(:,currOut);
            end
            ModelRuns = Sample.Post(:,StartInd:EndInd);
            postRuns =  Sample.PostPred(:,StartInd:EndInd);
            
            Model_at_PE_plus_noise = pointEstimatePredCurr(1,StartInd:EndInd) + postRuns - ModelRuns;
            
            % plot model output at sample
            %  postPlot = uq_violinplot(StartInd:EndInd,postRuns, 'FaceColor', postColor);
            postPlot = uq_violinplot(StartInd:EndInd,postRuns, 'FaceColor', postColor,...
                'FaceAlpha',0.0,    'EdgeColor',postColor,'LineWidth',3);
            hold on
            %                 Model_at_PE_plusPlot = uq_violinplot(StartInd:EndInd,...
            %                 Model_at_PE_plus_noise,'FaceColor',postModPEplusColor,'FaceAlpha',0.5);
            Model_at_PE_plusPlot = uq_violinplot(StartInd:EndInd,...
                Model_at_PE_plus_noise,'FaceColor',postModPEplusColor,...
                'FaceAlpha',0.0,    'EdgeColor',postModPEplusColor,'LineWidth',3,'LineStyle','--');
            
            % Update legend target
            legendObj = [postPlot(1) Model_at_PE_plusPlot(1)];
            legendName = {'posterior predictive',...
                ['(model at ',lower(pointEstimatePred{1}(1).Type),') +noise']};
            
        else
            error('mod_UQLab_plotSeriesPred: Plot option "post&Model(PE)+noise:edges" can only be used if point estimate is known');
        end % else to  if    if pointEstimatePredFlag
    otherwise
        fprintf('mod_UQLab_plotSeriesPred:\n')
        fprintf('Addmissible values for  plotType are: ''prior&Post'', ''prior'',''post'',\n')'
        fprintf(' ''post&Model'',  ''Model(PE)+noise'', ''post&Model(PE)+noise and ''.');
        error(' ''post&Model(PE)+noise:edges''.\n The used value ``%s`` is not valid,\n\n',plotType)
        
        
        
end




% Use scatter plot for observed data

DataY = Data.y;
xDummy = StartInd:EndInd;
for ii = xDummy
    yCurr = DataY(:,ii);
    for jj = 1:numel(yCurr)
        dataPlot = scatter(ii, yCurr(jj), 100, 'gx');
    end
end

% Update legend target
legendObj = [legendObj dataPlot(1)];
legendName{end+1} = 'data';

% Add point estimate
if pointEstimatePredFlag
    % assign model runs to current Data based on MOMap
    % define plot color order
    plotColors = uq_colorOrder(length(pointEstimatePred)+1);
    marker = '+';
    for pp = 1:length(pointEstimatePred)
        
        pointEstimatePredCurr = zeros(size(pointEstimatePred{pp}(1).Out,1),size(Data.y,2));
        % loop over data points
        for ii = 1:size(pointEstimatePredCurr,2)
            currModel = Data.MOMap(1,ii);
            currOut = Data.MOMap(2,ii);
            pointEstimatePredCurr(:,ii) = pointEstimatePred{pp}(currModel).Out(:,currOut);
        end
        
        %         for jj = 1:size(pointEstimatePredCurr,1)
        %             pointEstimatePlot = scatter(xDummy, pointEstimatePredCurr(jj,xDummy), 100, Marker,'MarkerEdgeColor',plotColors(1+pp,:));
        %
        %         end
        %
        
        % pointEstimatePlot = scatter(xDummy, pointEstimatePredCurr(:,xDummy), 100, 'r+');
        for jj = 1:size(pointEstimatePredCurr,1)
            pointEstimatePlot = scatter(xDummy, pointEstimatePredCurr(jj,xDummy), 100, marker,'MarkerEdgeColor',plotColors(1+pp,:));
        end
        % Update legend target
        legendObj = [legendObj,pointEstimatePlot(1)];
        legendName{end+1} = ['model at ',lower(pointEstimatePred{pp}(1).Type)];
        marker = '<';
    end
end
hold off

% Add legend
uq_legend(legendObj,legendName)

% Add axes labels
ylabel('$\mathcal{Y}$', 'Interpreter', 'LaTeX')
xlabel('$\mathrm{Data\,index\,}$(-)', 'Interpreter', 'LaTeX')

% Set xtick labels (maximum 10)
nTicks = 10;
if length(xDummy) > nTicks
    xDummy = ceil(linspace(1, length(xDummy), nTicks))-1+StartInd;
end
set(gca, 'XTick', xDummy)

% If ticklabel interpreter can be set, update labels to latex
labels = cell(length(xDummy),1);
if isprop(gca,'TickLabelInterpreter')
    for ii = 1:length(xDummy)
        labels{ii} = sprintf('$\\mathrm{y_{%u}}$',xDummy(ii));
    end
else
    for ii = 1:length(xDummy)
        labels{ii} = sprintf('y%d',xDummy(ii));
    end
end
set(gca, 'XTickLabel', labels)

end
