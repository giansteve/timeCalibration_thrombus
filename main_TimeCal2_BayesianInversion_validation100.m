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
root_destination = pwd;
addpath(root_destination)

% initialize UQlab
% uqlab

%% Load RawData
load('M:\IFM\User\melito\Phd\Projects\TimeCalibration_storageNoGitHub_saveFiles\Plot_AliModel_Calibration2_valid100\RawOutput_06June22_100sim.mat')

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
folderPath = 'TimeCal2_MRIdata';
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel2_Validation_final100\\%s',folderPath);
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_plotFittingMRI(fittedData,human_thr,'MRI_fitting')
cd(root_destination)


%% Load output: H/S and L/S
phic_HS_threshold = outToPCE.H_S;
phic_LS_threshold = outToPCE.L_S;

%% Define Input variables names
var_names = {'${D}_{\mathrm{c}}$',...
    '$\langle{\dot{\gamma}}\langle_\mathrm{t} $'};

%% Statistics on the input
fprintf('Get statistics INPUT ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel2_Validation_final100\\StatsInput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
figure('Visible','off')
plotmatrix(exp_design(:,[1 2]))
GM_printBMP(600,600,'InputStats')
GM_printEPS(600,600,'InputStats')
close all
cd(root_destination)

%% Statistics on the output
fprintf('Get statistics OUTPUT ... \n')
cd(root_destination)
try
    dest_plot = sprintf('Plot_AliModel2_Validation_final100\\StatsOutput');
    cd(dest_plot)
catch
    mkdir(dest_plot)
    cd(dest_plot)
end
timeCal_statisticsOutput(phic_HS_threshold,phic_LS_threshold,human_thr);
cd(root_destination)

%%
MRI_time = [1 11 21 31 41 51 61];
figure
% boxplot(phic_LS_threshold(MRI_time,:)');
hold on
% plot(human_thr.L_S(:,1))
plot(human_thr.time,human_thr.L_S(:,1),'r.','MarkerSize',10)
errorbar(human_thr.time,human_thr.L_S(:,1),human_thr.L_S(:,2),human_thr.L_S(:,2),'LineStyle','none','Marker','none','Color','r')


%%
figure
histogram(phic_LS_threshold(MRI_time(2),:))
% hold on
% histogram(phic_LS_threshold(3,:))
% histogram(phic_LS_threshold(4,:))
% histogram(phic_LS_threshold(5,:))
% histogram(phic_LS_threshold(6,:))
% histogram(phic_LS_threshold(7,:))

%%
[M,I] = (max(phic_LS_threshold(end,:)));
phic_HS_threshold_rmout = phic_HS_threshold;
phic_LS_threshold_rmout = phic_LS_threshold;
phic_HS_threshold_rmout(:,I) = [];
phic_LS_threshold_rmout(:,I) = [];
figure
subplot(211)
plot(linspace(0,1,61),phic_HS_threshold_rmout,'Color',[0.5 .5 .5])
hold on
plot(linspace(0,1,7),human_thr.H_S(:,1),'k.','MarkerSize',10)
errorbar(linspace(0,1,7),human_thr.H_S(:,1),human_thr.H_S(:,2),human_thr.H_S(:,2),'LineStyle','none','Marker','none','Color','k')
ylabel('$H/S$')
subplot(212)
plot(linspace(0,1,61),phic_LS_threshold_rmout,'Color',[0.5 .5 .5])
hold on
plot(linspace(0,1,7),human_thr.L_S(:,1),'k.','MarkerSize',10)
errorbar(linspace(0,1,7),human_thr.L_S(:,1),human_thr.L_S(:,2),human_thr.L_S(:,2),'LineStyle','none','Marker','none','Color','k')
ylabel('$L/S$')
xlabel('$t^*$ (-)')
GM_printBMP(500,300,'BiforcationOutput')
%%
[r,c] = find(phic_LS_threshold(end,:)<4.3);
firstGroupLS = phic_LS_threshold(MRI_time,c);

ED_firstGroup = exp_design(c,:);

[r,c] = find(phic_LS_threshold(end,:)<8);
secondGroupLS = phic_LS_threshold(MRI_time,c);
[r,c] = find(phic_LS_threshold(end,:)>4.3);
secondGroupLS = secondGroupLS(:,c);

ED_secondGroup = exp_design(c,:);
%%
colorStr = {'r','b'};
figure
for subplotCount = 1:2
    for m = 1:M
        if subplotCount == 1
            Y = ED_firstGroup(:,m);
        else
            Y = ED_secondGroup(:,m);
        end
        subplot(1,2,m)
        % The width of a histogram element is computed by the Scott's rule
        w = 3.49*std(Y)*numel(Y)^(-1/3);  % Width of a histogram element
        nBins = max(ceil(range(Y)/w),1);     % Number of histograms
        [hY,hX] = hist(Y,nBins);
        plot(hX,smooth(smooth(hY)),colorStr{subplotCount})
        hold on
    end
end
%%
figure
plotmatrix(ED_firstGroup)
figure
plotmatrix(ED_secondGroup)




















