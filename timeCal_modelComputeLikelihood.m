function [output] = timeCal_modelComputeLikelihood(X,P)
% function for the computation of the likelihood of agreement of fitting
% given the new set of input data X

% set script variables
% time_array = [1 11 21 31 41 51 61]; % indices for model instances
% fittedData = P.MRI_fittedData;
% PCE_H_S = P.PCE_H_S;
% PCE_L_S = P.PCE_L_S;
PCE_likelihood = P.PCE_likelihood;
% N = 1;

%% compute forward model
% H/S
% H_S_pce_eval = uq_evalModel(PCE_H_S,X);
% % L/S
% L_S_pce_eval = uq_evalModel(PCE_L_S,X);
% likelihood
output = uq_evalModel(PCE_likelihood,X);
% output
% Compute data difference
% diff_data.diff_H_S = H_S_pce_eval' - fittedData(time_array,2);
% diff_data.diff_L_S = L_S_pce_eval' - fittedData(time_array,3);
% %% compute likelihood
% diff_data.likelihoodH_S = ((var(diff_data.diff_H_S,[],1)).^(-N));
% diff_data.likelihoodL_S = ((var(diff_data.diff_L_S,[],1)).^(-N));
% diff_data.likelihoodH_S = diff_data.likelihoodH_S./sum(diff_data.likelihoodH_S);
% diff_data.likelihoodL_S = diff_data.likelihoodL_S./sum(diff_data.likelihoodL_S);
% % final likelihood
% weights_likelihood_HS_LS = [3 7]./10;
% likelihood_new = (weights_likelihood_HS_LS(1)./diff_data.likelihoodH_S.^(-1/N) + weights_likelihood_HS_LS(2)./diff_data.likelihoodL_S.^(-1/N) ).^N;
% likelihood_new = transpose(likelihood_new./sum(likelihood_new));
% output = likelihood_new;


end