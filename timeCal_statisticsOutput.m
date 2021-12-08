function [stats_prjct,output_matrix,diff_signal] = timeCal_statisticsOutput(output_matrix,exp_design,fittedData)
%MORPH_STATISTICSOUTPUT It makes the environment to create the statistical
%analysis of the output produced for the "TimeCalibration" project
% INPUT:
%       - output_matrix: it is the structure that includes all the
%           considered output for the SA

Ns = size(exp_design,1);
warning('off','all')
%% H/S thrombus
[stats_prjct.H_S] = timeCal_ConfInt(output_matrix.H_S_model,Ns);
% iOpts.Inference.Data = output_matrix.H_S_model';
% iOpts.Copula.Type = 'Independent';
% stats_prjct.H_S.pdf = uq_createInput(iOpts);

%% L/S thrombus
[stats_prjct.L_S] = timeCal_ConfInt(output_matrix.L_S_model,Ns);
% iOpts.Inference.Data = output_matrix.L_S_model';
% iOpts.Copula.Type = 'Independent';
% stats_prjct.L_S.pdf = uq_createInput(iOpts);
%% Surface area
% [stats_prjct.SA] = timeCal_ConfInt(output_matrix.SA_model,Ns);
% iOpts.Inference.Data = output_matrix.SA_model';
% iOpts.Copula.Type = 'Independent';
% stats_prjct.SA.pdf = uq_createInput(iOpts);
%% Volume
% [stats_prjct.Vol] = timeCal_ConfInt(output_matrix.VOL_model,Ns);
% iOpts.Inference.Data = output_matrix.VOL_model';
% iOpts.Copula.Type = 'Independent';
% stats_prjct.Vol.pdf = uq_createInput(iOpts);
%% Plot
diff_signal = timeCal_statsPlot(stats_prjct,output_matrix,fittedData);
timeCal_plotFinalResult(output_matrix,fittedData);
%% REMEMBER TO UNCOMMENT
warning('on','all')
close all
end
















