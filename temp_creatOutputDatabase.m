% %% Output
% folderPath = 'TimeCal_Out';
% % Collect output from the solution files
% fprintf('Reading the output files ... \n')
% [OUTPUT,crushed_sim_idx] = timeCal_read_output(Ns,folderPath);
% 
% % delete crushed simulations and relative exp_design
% OUTPUT(crushed_sim_idx) = [];
% exp_design(crushed_sim_idx,:) = [];
% 
% % Transform the structure into sth readable
% fprintf('Transforming output ... \n')
% thresholdThr = 0.50;
% [outToPCE] = timeCal_getOutputReadable(OUTPUT,thresholdThr);
% 
% % safe
% fprintf('Saving the workspace ... \n')
% save('RawOutput_23Feb21_1304.mat','-v7.3')

%% ===================================================================
% load('RawOutput_23Feb21_1304.mat')
OUT_time = outToPCE.time(:,1);
OUT_HS = outToPCE.H_S(:);
OUT_LS = outToPCE.L_S(:);
save('TimeCal_OUTPUT.mat')
%% ===================================================================