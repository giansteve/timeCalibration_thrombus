function [RA_set] = timeCal_RAmain(output,model,Ns)
% set and compute the reliability analysis on 

%% Define INPUT module
% [~,idx_cache] = sort(SA_strct.total);
% nonConst = idx_cache(end-3:end);
nonConst = 1:10;
Input_RA = timeCal_createInput(Ns,nonConst,0,[]);
% relopts.Input = INPUT_RA;

%% set model for computing RA
modelOpts.mFile = 'timeCal_RAmodel';
modelOpts.Parameters.mod = model;
modelOpts.Parameters.threshold = output.q25_data;
RA_model = uq_createModel(modelOpts);
% RA_modelEval = uq_evalModel(RA_model,exp_design_RA);

%% MCS
reloptsRA.Type = 'uq_reliability';
reloptsRA.Method = 'MCS';
reloptsRA.Display = 'quiet';
reloptsRA.Model = RA_model;
reloptsRA.Simulation.MaxSampleSize = Ns;
RA_set = uq_createAnalysis(reloptsRA);


% %% AKMCS
% relopts.Type = 'uq_reliability';
% relopts.Method = 'AKMCS';
%  
% 
% relopts.Model = RA_model;
% 
% % get corrected output exp_design
% % output_RA = zeros(size(exp_design,1),1);
% % output_RA(output.indicator <= output.q25_data) = 1;
% 
% % set threshold
% relopts.LimitState.CompOp = '<';
% relopts.LimitState.Threshold = 0;
% 
% % set method
% % relopts.Simulation.BatchSize = 10000; % # samples evaluated at once
% 
% % set AKMCS
% relopts.AKMCS.MetaModel = 'kriging';
% relopts.AKMCS.PCK.Mode = 'sequential';
% % relopts.AKMCS.PCK.Input = relopts.Input;
% % relopts.AKMCS.PCK.ExpDesign.X = exp_design;
% % relopts.AKMCS.PCK.ExpDesign.Y = output_RA;
% % relopts.AKMCS.PCK.PCE.Degree = 5;
% relopts.AKMCS.PCK.Kriging.Corr.Family = 'Gaussian';
% relopts.AKMCS.Convergence = 'stopU';
% relopts.AKMCS.MaxAddedED = 100;
% relopts.AKMCS.IExpDesign.N = length(exp_design_RA);
% relopts.AKMCS.IExpDesign.Sampling = 'lhs';
% % reducing exp_design to nonConst RVs
% % relopts.AKMCS.IExpDesign.X = exp_design_RA(:,relopts.Input.nonConst);
% % relopts.AKMCS.IExpDesign.X = exp_design_RA;
% % relopts.AKMCS.IExpDesign.G = RA_modelEval;
% 
% % storage evaluation
% % relopts.SaveEvaluation = 1;
% 
% 
% 
% RA_set = uq_createAnalysis(relopts);






end % function