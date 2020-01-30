function [ data ] = LoadData(varargin)
%LoadData Extracts data from the APS .dat files
%   Gets the data from the APS scan files. 
%   First arg is file location from current.
%   Second arg is the base file name.
%   Third arg is the scan number.

%Gets path name
pathName = pwd;

switch nargin
    case 3
        location = varargin{1};
        Name = varargin{2};
        scan = varargin{3};
        %Checks for \ and . (adds them if missing)
        if ~strcmp(location(1),'\')
            location = strcat('\',location);
        else
           location = varargin{1};
        end
        if ~strcmp(location(end),'\')
            location = strcat(location,'\');
        end
        if ~strcmp(Name(end),'.')
            Name = strcat(Name,'.');
        end
        
        fileName = strcat(location,Name,num2str(sprintf('%04d',scan)),'.dat');
        data.Name = strcat(fileName,num2str(sprintf('%04d',scan)),'.dat');
        data.Path = strcat(pathName,location);
        
    case 1
        
        fileName = varargin{1};
end

%Creates fullFileName if not equal
fullFileName = strcat(pathName,fileName);

%Saves the Name, fullFileName and pathName into the dat
data.Full = fullFileName;

%Opens the file to identifier fid
fid = fopen(data.Full);

%Looks for Scan number in file
allFile = textscan(fid,'%s','Delimiter','\n'); %Saves all of the file to a giant string
row = find(~cellfun('isempty',strfind(allFile{1},'Scan number'))); %Gets the number of the row from the cell
pos = strfind(allFile{1}{row},' '); %Gets the positions of the qotations
data.scan = round(str2double(allFile{1}{row}(pos(4)+1:end))); %Gets the Bfield and stores it

%Looks for Temp in file
row = find(~cellfun('isempty',strfind(allFile{1},'4idc1:LS340:TC2:SP1'))); %Gets the number of the row from the cell
if isempty(row)
    print('Warning! Temperature data not found for file: ', fileName{n});
end
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
data.T = round(str2double(allFile{1}{row}(pos(1)+1:pos(2)-1))); %Gets the Bfield and stores it

%Looks for Bfield in file
row = find(~cellfun('isempty',strfind(allFile{1},'4idc2:AMI430:Field'))); %Gets the number of the row from the cell
if isempty(row)
    print('Warning! Magnetic Field data not found for file: ', fileName{n});
end
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
data.B = round(str2double(allFile{1}{row}(pos(1)+1:pos(2)-1)),4); %Gets the Bfield and stores it

%Looks for Bfield in file
row = find(~cellfun('isempty',strfind(allFile{1},'4idc2:AMI430:Current'))); %Gets the number of the row from the cell
if isempty(row)
    print('Warning! Magnetic Field data not found for file: ', fileName{n});
end
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
data.Bcurrent = round(str2double(allFile{1}{row}(pos(1)+1:pos(2)-1))/10,4); %Gets the Bfield and stores it

%Looks for zPos in file
try
row = find(~cellfun('isempty',strfind(allFile{1},'4idc1:m13.RBV'))); %Gets the number of the row from the cell
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
data.zPos = str2double(allFile{1}{row}(pos(1)+1:pos(2)-1)); %Gets the Bfield and stores it
catch
    %warning('Could not find Z Position.');
end

%Looks for theta in file
try
row = find(~cellfun('isempty',strfind(allFile{1},'4idc1:m19.RBV'))); %Gets the number of the row from the cell
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
data.theta = str2double(allFile{1}{row}(pos(1)+1:pos(2)-1)); %Gets the Bfield and stores it
catch
    %warning('Could not find theta Position.');
end

%Gets i0 sensitivity
row = find(~cellfun('isempty',strfind(allFile{1},'4idc1:A1sens_num.VAL'))); %Gets the number of the row from the cell
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
num = str2double(allFile{1}{row}(pos(1)+1:pos(2)-1)); %Gets the sensitivity number and stores it
row = find(~cellfun('isempty',strfind(allFile{1},'4idc1:A1sens_unit.VAL'))); %Gets the number of the row from the cell
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
unit = allFile{1}{row}(pos(1)+1:pos(2)-1); %Gets the unit and stores it
switch unit
    case 'pA/V'
        data.sensitivity.i0 = num*10^-12;
    case 'nA/V'
        data.sensitivity.i0 = num*10^-9;
    case 'uA/V'
        data.sensitivity.i0 = num*10^-6;        
end

%Gets TEY sensitivity
row = find(~cellfun('isempty',strfind(allFile{1},'4idc1:A6sens_num.VAL'))); %Gets the number of the row from the cell
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
num = str2double(allFile{1}{row}(pos(1)+1:pos(2)-1)); %Gets the sensitivity number and stores it
row = find(~cellfun('isempty',strfind(allFile{1},'4idc1:A6sens_unit.VAL'))); %Gets the number of the row from the cell
pos = strfind(allFile{1}{row},'"'); %Gets the positions of the qotations
unit = allFile{1}{row}(pos(1)+1:pos(2)-1); %Gets the unit and stores it
switch unit
    case 'pA/V'
        data.sensitivity.tey = num*10^-12;
    case 'nA/V'
        data.sensitivity.tey = num*10^-9;
    case 'uA/V'
        data.sensitivity.tey = num*10^-6;        
end


fseek(fid,0,'bof'); %Resets position to beginning of file.
%%%%%%%%%%%%%%%%%



%Pulls the data from the file
data.raw{1} = [];
while isempty(data.raw{1})
    data.raw = textscan(fid, '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f %*[^\n]','headerlines',1);
end;
%Closes the file
fclose(fid);

end

