function [] = model_funct_TimeCal(X)
%MODEL_FUNCT here we copy the folder for the n-th simulation, identify and
%write on the input file, start simulation and get the output.

for sim = 1:size(X,1)
    
    %     time = tic;
    
    new_name = sprintf('TC_X%.4d/',sim); % new name of folder
    copyfile('0_0_1_xuPimpleFoamThrombusBFSSAPeriod/',new_name); % make the new simulation folder
    
    cd(new_name); % directory to new folder
    
    %% make jobscript file for each sim
    fileID = fopen("jobscript",'r+');
    dummyline = fgetl(fileID); % get the 1st line
    counterj = 1;
    
    while ~strcmpi(dummyline(1:16),'simulation_path=')
        dummyline = fgetl(fileID);
        counterj = counterj + 1 ;
        while length(dummyline) < 3 || isempty(dummyline)
            dummyline = fgetl(fileID);
            counterj = counterj + 1 ;
        end
    end
    oldline = dummyline;
    fileID = fopen("jobscript",'r+');
    for line = 1:counterj-1
        dummyline = fgetl(fileID);
    end
    oldvalue =  "/home/azerila/OpenFOAM/azerila-5.x/run/callibration/BFS_Newtonian";
    str2 =convertCharsToStrings(new_name);
    newvalue = fullfile(oldvalue,str2);
    %     newvalue2 = sprintf('%s/%s%s\n',oldvalue,str2,' ');
    %     newvalue =  " /home/azerila/OpenFOAM/azerila-5.x/run/callibration/BFS_Newtonian ",new_name
    Y = strrep(oldline,oldvalue,newvalue);
    fwrite(fileID,Y); % write on file
    fclose(fileID);
    
    
    %% Input file
    inputFile_Dc = 'xutransportProperties';
    input_Dc(inputFile_Dc,X(sim,1));
    
    inputFile_Kc = 'xutransportProperties';
    input_Kc(inputFile_Kc,X(sim,2));
    
    inputFile_K_BP = 'xutransportProperties';
    input_K_BP(inputFile_K_BP,X(sim,3));
    
    inputFile_Ct = 'xutransportProperties';
    input_Ct(inputFile_Ct,X(sim,4));
    
    inputFile_BPt = 'xutransportProperties';
    input_BPt(inputFile_BPt,X(sim,5));
    
    inputFile_RTt = 'xutransportProperties';
    input_RTt(inputFile_RTt,X(sim,6));
    
    inputFile_KCwall = 'xutransportProperties';
    input_KCwall(inputFile_KCwall,X(sim,7));
    
    inputFile_C_AP = 'xutransportProperties';
    input_C_AP(inputFile_C_AP,X(sim,8));
    
    inputFile_BPbt = 'xutransportProperties';
    input_BPbt(inputFile_BPbt,X(sim,9));
    
    inputFile_GammaDot = 'switches';
    input_GammaDot(inputFile_GammaDot,X(sim,10));
    
    %% simulation starts in 3... 2... 1...
    
    %     system('./Allclean');
    %     system('./Allrun');
    
    %% Get the output
    %     outputFileName = 'postProcessing';
    %
    %     temp_result = read_output(outputFileName);
    %
    %     if size(temp_result,1) ~= 100
    %         result(sim,:) = repmat(-1,100,1);
    %     else
    %         result(sim,:) = temp_result;
    %     end
    
    cd ../
    
    %     time = toc(time);
    %     txt = sprintf('Simulation time is: %3f',time);
    %     fprintf(txt)
    %
    %     pause(0.01)
    %
    % plot
    %     figure(1)
    %     plot(linspace(1,100,100),result(sim,:))
    %     hold on
    %     grid on
    %     title('Results')
    %     xlabel('Time [s]')
    %     ylabel('Vol average of thrombus')
    
    
end



end

