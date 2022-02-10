function [OUTPUT,crushed_sim_idx] = timeCal_read_output(Ns,folderPath)
% Read and report the output of the morphological_aorta project.

%% destination
dest_folder = sprintf('C:\\Users\\gm20m18\\Desktop\\%s',folderPath);
cd(dest_folder)
%% reading data
% foldersSim = dir(pwd);                      % get the current folder
% foldersSim([foldersSim.isdir] == 0) = [];   % delete non folders
% foldersSim(1:2) = [];                       % delete going back folders
% clear foldersSim

crush_idx = 1;
thresholdThr = 0.98;
crushed_sim_idx = [];
% date_begin = datetime('now');

for nsim = 1:Ns
    text_screen_sim = sprintf(' - simulation: %d/%d: ',nsim,Ns);
%     elaps_time(nsim) = tic;
    %     cd(sprintf('%s\\%s\\postProcessing',pwd,foldersSim(nsim).name));
    
    try % is the simulation XXXX done?

        if nsim < 10 && nsim >= 1
            cd(sprintf('%s\\TC_X000%d\\postProcessing',pwd,nsim));
        elseif nsim <= 99 && nsim >= 10
            cd(sprintf('%s\\TC_X00%d\\postProcessing',pwd,nsim));
        elseif nsim <= 999 && nsim >= 100
            cd(sprintf('%s\\TC_X0%d\\postProcessing',pwd,nsim));
        else
            cd(sprintf('%s\\TC_X%d\\postProcessing',pwd,nsim));
        end
        
    catch % simulation not performed
        warning('off','all')
        txt_sim = sprintf('\t NO :( - \n');
        fprintf(txt_sim)
        crushed_sim_idx(crush_idx) = nsim;
        crush_idx = crush_idx + 1;
        cd(dest_folder)
        continue
    end
    
    
    try % in case of crushed sim
        fprintf(text_screen_sim)
        
        %% H/S
        % Coordinates
        cd('Probes_faceCentres_STEP_H\0')
        fileID_phiC = fopen('phi_c');
        filePhiC_scan = textscan(fileID_phiC,'%s','Delimiter','\n');
        probeCoord = find(cell2mat(strfind(filePhiC_scan{1,1},'Probe')) == 3);
        OUTPUT(nsim).H_S.Tot_pts = size(probeCoord,1);
        for probe_point = 1:OUTPUT(nsim).H_S.Tot_pts
            probCoord_temp{probe_point,1} = strsplit(filePhiC_scan{1,1}{probe_point});
            OUTPUT(nsim).H_S.Coordinate(probe_point,:) = [str2doubles(probCoord_temp{probe_point,1}{1,4}(2:end)) str2doubles(probCoord_temp{probe_point,1}{1,5})];
        end
        OUTPUT(nsim).H_S.Length = max(OUTPUT(nsim).H_S.Coordinate) - min(OUTPUT(nsim).H_S.Coordinate);
        fclose(fileID_phiC);
        % VFT reading
        fileID_phiC = fopen('phi_c');
        H_S_phiC_temp = textscan(fileID_phiC,'%s','HeaderLines',OUTPUT(nsim).H_S.Tot_pts+2,'Delimiter','\t');
        for time_idx = 1:numel(H_S_phiC_temp{1,1})
            OUTPUT(nsim).H_S.phi_c(time_idx,:) = str2doubles(strsplit(H_S_phiC_temp{1,1}{time_idx}));
            logic_temp = OUTPUT(nsim).H_S.phi_c(time_idx,2:end) >= thresholdThr;
            if sum(logic_temp) ~= 0
                OUTPUT(nsim).H_S.output(time_idx,1) = (max(find(logic_temp)) / OUTPUT(nsim).H_S.Tot_pts * OUTPUT(nsim).H_S.Length(2)) / OUTPUT(nsim).H_S.Length(2);
            else
                OUTPUT(nsim).H_S.output(time_idx,1) = 0;
            end
        end
        fclose(fileID_phiC);
        OUTPUT(nsim).time = OUTPUT(nsim).H_S.phi_c(:,1);
        if size(OUTPUT(nsim).time,1) ~= 80 % sim didnt perform the whole time?
            warning('off','all')
            txt_sim = sprintf('\t NO :( - \n');
            fprintf(txt_sim)
            crushed_sim_idx(crush_idx) = nsim;
            crush_idx = crush_idx + 1;
            cd(dest_folder)
            continue
        end
        cd ..\..
        
        %% L/S
        % Coordinates
        cd('Probes_faceCentres_STEP_L\0')
        fileID_phiC = fopen('phi_c');
        filePhiC_scan = textscan(fileID_phiC,'%s','Delimiter','\n');
        probeCoord = find(cell2mat(strfind(filePhiC_scan{1,1},'Probe')) == 3);
        OUTPUT(nsim).L_S.Tot_pts = size(probeCoord,1);
        for probe_point = 1:OUTPUT(nsim).L_S.Tot_pts
            probCoord_temp{probe_point,1} = strsplit(filePhiC_scan{1,1}{probe_point});
            OUTPUT(nsim).L_S.Coordinate(probe_point,:) = [str2doubles(probCoord_temp{probe_point,1}{1,4}(2:end)) str2doubles(probCoord_temp{probe_point,1}{1,5})];
        end
        OUTPUT(nsim).L_S.Length = max(OUTPUT(nsim).L_S.Coordinate) - min(OUTPUT(nsim).L_S.Coordinate);
        fclose(fileID_phiC);
        % VFT reading
        fileID_phiC = fopen('phi_c');
        L_S_phiC_temp = textscan(fileID_phiC,'%s','HeaderLines',OUTPUT(nsim).L_S.Tot_pts+2,'Delimiter','\t');
        for time_idx = 1:numel(L_S_phiC_temp{1,1})
            OUTPUT(nsim).L_S.phi_c(time_idx,:) = str2doubles(strsplit(L_S_phiC_temp{1,1}{time_idx}));
            logic_temp = OUTPUT(nsim).L_S.phi_c(time_idx,2:end) >= thresholdThr;
            if sum(logic_temp) ~= 0
                OUTPUT(nsim).L_S.output(time_idx,1) = (max(find(logic_temp)) / OUTPUT(nsim).L_S.Tot_pts * OUTPUT(nsim).L_S.Length(1)) / OUTPUT(nsim).H_S.Length(2);
            else
                OUTPUT(nsim).L_S.output(time_idx,1) = 0;
            end
        end
        fclose(fileID_phiC);
        cd ..\..
        
        %% Surface area
        %         cd('swakExpression_SurefaceArea\0')
        %         fileID_SA = fopen('SurefaceArea');
        %         fileSA_scan = textscan(fileID_SA,'%s','Delimiter','\n','HeaderLines',1);
        %         for time_idx = 1:numel(fileSA_scan{1,1})
        %             OUTPUT(nsim).SA.output(time_idx,:) = str2doubles(strsplit(fileSA_scan{1,1}{time_idx}));
        %         end
        %         OUTPUT(nsim).SA.output(:,1) = [];
        %         fclose(fileID_SA);
        %         cd ..\..
        
        %% Volume
        %         cd('volumeAverage_phi_c_vol_average\0.01')
        %         fileID_phiC = fopen('phi_c');
        %         fileVol_scan = textscan(fileID_phiC,'%s','Delimiter','\n','HeaderLines',1);
        %         for time_idx = 1:numel(fileVol_scan{1,1})
        %             OUTPUT(nsim).Vol.output(time_idx,:) = str2doubles(strsplit(fileVol_scan{1,1}{time_idx}));
        %         end
        %         OUTPUT(nsim).Vol.output(:,1) = [];
        %         fclose(fileID_phiC);
        %         cd ..\..
        
        cd ..\.. % go back to sims folder
        warning ('off','all')
        txt_sim = sprintf('\t OK :) \n');
        fprintf(txt_sim)
        warning ('on','all')
        
    catch % sth went wrong during simulation
        txt_sim = sprintf('\t NO :( - \n');
        warning('off','all')
        fprintf(txt_sim)
        crushed_sim_idx(crush_idx) = nsim;
        crush_idx = crush_idx + 1;
        cd(dest_folder)
        continue
    end % try loop for crushed sims
    fclose all;
end % for loop on simulations


cd ..\..\..

end % function


