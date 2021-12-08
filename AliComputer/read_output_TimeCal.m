function [result] = read_output_TimeCal(outputFileName)
%READ_OUTPUT Summary of this function goes here
%   Detailed explanation goes here

cd('./0.2');


fileID = fopen(outputFileName);

for ll = 1:19
    dummyline = fgetl(fileID);
end

temp_res = fgetl(fileID);
result = str2num(temp_res(24:end-1));
fclose(fileID)

cd ../..

end

