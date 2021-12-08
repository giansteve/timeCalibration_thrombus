function [] = timeCal_statisticsOutput(OUT_HS,OUT_LS,OUT_time,fittedData)


% Ns = size(exp_design,1);
% warning('off','all')

% No outliers detection and deletion
%% H/S thrombus
% [stats_prjct.H_S] = timeCal_ConfInt(output_matrix.H_S_model,Ns);

% %% L/S thrombus
% [stats_prjct.L_S] = timeCal_ConfInt(output_matrix.L_S_model,Ns);

%% Plot
timeCal_statsPlot(OUT_HS,OUT_LS,OUT_time);
timeCal_plotFinalResult(OUT_HS,OUT_LS,OUT_time,fittedData);

% warning('on','all')
close all
end
















