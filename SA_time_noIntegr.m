function [S_cum,ST_cum] = SA_time_noIntegr(uq_pce_struct,M)
% This function is able to get the indices and coefficient of the
% polynomial chaos expansion that has been derived from the UQlab toolbox
% and get the related Sensitivity Analysis indices from Sobol'. The
% particular formulation is related to time series. 
% THIS VERSION TRIES NOT TO INTEGRATE THE SOBOL INDICES IN TIME

% M = length(uq_pce_struct.Internal.Input.Marginals);
% Computing elements for Sobol's indices

% pre-allocating
S_cum = zeros(length(uq_pce_struct),M);
ST_cum = zeros(length(uq_pce_struct),M);
% v_i = S_cum;
% v_tot = S_cum;
% v_set = zeros(length(uq_pce_struct),1);

m_seq = 1:M;
for t = 1:length(uq_pce_struct)
    for m = 1:M
        alpha_m_vec = uq_pce_struct(t).Basis.Indices(:,m_seq(m_seq~=m));
        [A_i_row,~] = find(all(alpha_m_vec == 0,2));
        [A_i_tot,~] = find(uq_pce_struct(t).Basis.Indices(:,m) > 0);
        % get variance given by the input m - no interaction terms
        S_cum(t,m) = sum(uq_pce_struct(t).Coefficients(A_i_row(2:end)).^2);
        % get variance given by the input m - yes interaction terms
        ST_cum(t,m) = sum(uq_pce_struct(t).Coefficients(A_i_tot).^2);
    end
    % get total variance at time step t
%     v_set(t) = sum(uq_pce_struct(t).Coefficients(2:end).^2);
end

% Computing Sobol'indices
% Denom = v_set;
% for m = 1:M
%     S_cum(:,m) = (v_i(:,m))./Denom';
%     ST_cum(:,m) = (v_tot(:,m))./Denom';
% end
% S_cum(Denom == 0,:) = 0;
% ST_cum(Denom == 0,:) = 0;




end