function [SA] = timeCal_SA(PCE,M,var_names)
% Perform sensitivity analysis on the TimeCalibration project


% H/S
[SA.H_S.main,SA.H_S.total] = SA_time(PCE.H_S.time,M);
% L/S
[SA.L_S.main,SA.L_S.total] = SA_time(PCE.L_S.time,M);
% SA
[SA.SA.main,SA.SA.total] = SA_time(PCE.SA.time,M);
% Vol
[SA.Vol.main,SA.Vol.total] = SA_time(PCE.VOL.time,M);

% diff Signal H/S                                         
[SA.diff_indicat.H_S.main,SA.diff_indicat.H_S.total,SA.diff_indicat.H_S.names] = SA_Sobols(PCE.diff_indicat.H_S,M,var_names,0);
% diff Signal L/S
[SA.diff_indicat.L_S.main,SA.diff_indicat.L_S.total,SA.diff_indicat.L_S.names] = SA_Sobols(PCE.diff_indicat.L_S,M,var_names,0);
% diff Signal SA
[SA.diff_indicat.SA.main,SA.diff_indicat.SA.total,SA.diff_indicat.SA.names] = SA_Sobols(PCE.diff_indicat.SA,M,var_names,0);
% diff Signal Vol
[SA.diff_indicat.Vol.main,SA.diff_indicat.Vol.total,SA.diff_indicat.Vol.names] = SA_Sobols(PCE.diff_indicat.Vol,M,var_names,0);

% % Indicator diff Signal H/S
% [SA.diff_indicat.H_S.main,SA.diff_indicat.H_S.total] = SA_time(PCE.diff_indicat.H_S,M);
% % Indicator diff Signal L/S
% [SA.diff_indicat.L_S.main,SA.diff_indicat.L_S.total] = SA_time(PCE.diff_indicat.L_S,M);
% % Indicator diff Signal SA
% [SA.diff_indicat.SA.main,SA.diff_indicat.SA.total] = SA_time(PCE.diff_indicat.SA,M);
% % Indicator diff Signal Vol
% [SA.diff_indicat.Vol.main,SA.diff_indicat.Vol.total] = SA_time(PCE.diff_indicat.Vol,M);
end