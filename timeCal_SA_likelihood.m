function [SA] = timeCal_SA_likelihood(PCE,M,var_names)
% Perform sensitivity analysis on the TimeCalibration project


% H/S
[SA.H_S.main,SA.H_S.total] = SA_time_noIntegr(PCE.H_S,M);
% L/S
[SA.L_S.main,SA.L_S.total] = SA_time_noIntegr(PCE.L_S,M);
% % SA
% [SA.SA.main,SA.SA.total] = SA_time_noIntegr(PCE.SA,M);
% % Vol
% [SA.Vol.main,SA.Vol.total] = SA_time_noIntegr(PCE.VOL,M);

% LOGLikelihood
% [SA.LOGlikeL.main,SA.LOGlikeL.total,SA.LOGlikeL.names_sort] = SA_Sobols(PCE.diff_indicat.LOGlikelihood,M,var_names,0);

% LOG L/S likelihood
[SA.LOGlikeLS.main,SA.LOGlikeLS.total,SA.LOGlikeLS.names_sort] = SA_Sobols(PCE.diff_indicat.LSlogLikelihood,M,var_names,0);

end