function [] = input_Dc(input_fileName,X)
%INPUT_DC

cd ./constant/

fileID = fopen(input_fileName,'r+');
dummyline = fgetl(fileID); % get the 1st line
counter = 1;

while ~strcmpi(dummyline(1:3),'D_c')
    dummyline = fgetl(fileID);
    counter = counter + 1 ;
    while length(dummyline) < 3 || isempty(dummyline)
        dummyline = fgetl(fileID);
        counter = counter + 1 ;
    end
end
fclose(fileID);
oldline = dummyline;
fileID = fopen(input_fileName,'r+');
for line = 1:counter-1
    dummyline = fgetl(fileID);
end

oldvalue = '1.00e-08';
newvalue = num2str(X,'%1.2e');
Y = strrep(oldline,oldvalue,newvalue);
fwrite(fileID,Y); % write on file
fclose(fileID);

cd ../

end

