function [S_main_sort,S_TOT_sort,names_sort] = SA_Sobols(uq_pce_struct,M,var_names,sort_var)
% This function is able to get the indices and coefficient of the
% polynomial chaos expansion that has been derived from the UQlab toolbox
% and get the related Sensitivity Analysis indices from Sobol'. The
% particular formulation is related to time series. Such application has
% been not yet implemented into the UQlab toolbox.


% Computing elements for Sobol's indices
if size(uq_pce_struct,2) == 1 % if N_out = 1 % added: .time
    % pre-allocating
    S_main = zeros(length(uq_pce_struct),M);
    S_TOT = zeros(length(uq_pce_struct),M);
    v_i = S_main;
    v_tot = S_main;
    % Computing
    m_seq = 1:M;
    for m = 1:M
        alpha_m_vec = uq_pce_struct.Basis.Indices(:,m_seq(m_seq~=m));
        [A_i_row,~] = find(all(alpha_m_vec == 0,2));
        [A_i_tot,~] = find(uq_pce_struct.Basis.Indices(:,m) > 0);
        % get variance given by the input m - no interaction terms
        v_i(m) = sum(uq_pce_struct.Coefficients(A_i_row(2:end)).^2);
        % get variance given by the input m - yes interaction terms
        v_tot(m) = sum(uq_pce_struct.Coefficients(A_i_tot).^2);
    end
    % get total variance at time step t
    v_set = sum(uq_pce_struct.Coefficients(2:end).^2);
else % for multiple Outputs
    % pre-allocating
    S_main = zeros(size(uq_pce_struct,2),M);
    S_TOT = zeros(size(uq_pce_struct,2),M);
    v_i = S_main;
    v_tot = S_main;
    % Computing
    for n_out = 1:size(uq_pce_struct,2)
        m_seq = 1:M;
        for m = 1:M
            alpha_m_vec = uq_pce_struct(n_out).Basis.Indices(:,m_seq(m_seq~=m));
            [A_i_row,~] = find(all(alpha_m_vec == 0,2));
            [A_i_tot,~] = find(uq_pce_struct(n_out).Basis.Indices(:,m) > 0);
            % get variance given by the input m - no interaction terms
            v_i(n_out,m) = sum(uq_pce_struct(n_out).Coefficients(A_i_row(2:end)).^2);
            % get variance given by the input m - yes interaction terms
            v_tot(n_out,m) = sum(uq_pce_struct(n_out).Coefficients(A_i_tot).^2);
        end
        % get total variance at time step t
        v_set(n_out) = sum(uq_pce_struct(n_out).Coefficients(2:end).^2);
    end
end

% Computing Sobol'indices
% Denom = cumsum(v_set);
for n_out = 1:size(uq_pce_struct,2)
    for m = 1:M
        S_main(n_out,m) = cumsum(v_i(n_out,m))./v_set(n_out);
        S_TOT(n_out,m) = cumsum(v_tot(n_out,m))./v_set(n_out);
    end
end
names_sort = var_names;
S_main_sort = S_main;
S_TOT_sort = S_TOT;

if sort_var ~= 0 && size(uq_pce_struct,2) == 1
    [S_main_sort,idx_SA_sort] = sort(S_main,'descend');
    S_TOT_sort = S_TOT(idx_SA_sort);
    
    % ordering the sobol indices and getting the interaction terms
    for idx = 1:length(idx_SA_sort)
        names_sort{1,idx} = var_names{1,idx_SA_sort(idx)};
    end
    
end




end