function [Grid_N,maxlevel,Coarsening,K_original] = reader_Homogenizer(Directory, InputFile, PermDir)

%% Read Input File
File = strcat(Directory,'/',InputFile);
fileID = fopen(File, 'r');
% Read lines from input file
matrix = textscan(fileID, '%s', 'Delimiter', '\n');
InputMatrix = matrix{1};
fclose(fileID);
% Remove lines which are commented (contain --)
Commented = startsWith(InputMatrix, '--');
InputMatrix(Commented) = {'--'}; % removing the string if it is commented.

%% Read SettingFile
SettingsFile = strcat(Directory, '/SimulatorSettings.txt');
fileID = fopen(SettingsFile, 'r');
matrix = textscan(fileID, '%s', 'Delimiter', '\n');
SettingsMatrix = matrix{1};
fclose(fileID);
% Remove lines which are commented (contain --)
Commented = startsWith(SettingsMatrix, '--');
SettingsMatrix(Commented) = {'--'};

%% Read GRID
temp = strfind(InputMatrix, 'SPECGRID');
index = find(~cellfun('isempty', temp));
Grid_N = [str2double(InputMatrix{index+1});...
    str2double(InputMatrix{index+2});
    str2double(InputMatrix{index+3})];

% Read coarsening ratios
temp = strfind(SettingsMatrix, 'ADM');
adm = find(~cellfun('isempty', temp));
if str2double(SettingsMatrix(adm + 1)) == 1
    temp = strfind(SettingsMatrix, 'LEVELS');
    x = find(~cellfun('isempty', temp));
    maxlevel = str2double(SettingsMatrix(x+1));
    temp = strfind(SettingsMatrix, 'COARSENING_RATIOS');
    x = find(~cellfun('isempty', temp));
    Coarsening(1) = str2double(SettingsMatrix(x+1));
    Coarsening(2) = str2double(SettingsMatrix(x+2));
    Coarsening(3) = str2double(SettingsMatrix(x+3));
end

%% Read Perm
perm = zeros(3, 1);
temp = strfind(InputMatrix, 'PERMX');
perm(1) = find(~cellfun('isempty', temp));
temp = strfind(InputMatrix, 'PERMY');
perm(2) = find(~cellfun('isempty', temp));
temp = strfind(InputMatrix, 'PERMZ');
perm(3) = find(~cellfun('isempty', temp));
for i=1:3
    if strcmp(InputMatrix(perm(i) - 1), 'INCLUDE')
        PermInclude(i) = 1;
        PermFile{i} = strcat(PermDir,'/',char(InputMatrix(perm(i) +1)));
    else
        PermInclude(i) = 0;
        Perm(i) = str2double(InputMatrix(perm(i) +1));
    end
end
K_original = ones(Grid_N(1)*Grid_N(2)*Grid_N(3), 3);
for i=1:3
    if PermInclude(i)
        % load the file in a vector
        if i==1
            fprintf('\n---> Reading permeability file #%d ...',i);
            field = load(PermFile{i});
            fprintf(' ---> Completed.\n');
        else
            % loading the permeability file only if different
            % than previous one
            if ~strcmp(PermFile{i},PermFile{i-1})
                fprintf('---> Reading permeability file #%d ...',i);
                field = load(PermFile{i});
                fprintf(' ---> Completed.\n');
            end
        end
        % reshape it to specified size
        if (Grid_N(1)~=field(1,1)) || (Grid_N(2)~=field(2,1)) || (Grid_N(3)~=field(3,1))
            warning('The grid cells mentioned in the permeability file #%d do not match with Reservoir grid cells.\n',i);
        end
        if (Grid_N(1)>field(1,1)) || (Grid_N(2)>field(2,1)) || (Grid_N(3)>field(3,1))
            error('The grid cells mentioned in the permeability file #%d are less than that of Reservoir grid cells. Check the input file.\n',i);
        end
        field1 = reshape(field(4:end,1),[field(1,1) field(2,1) field(3,1)]);
        % make it the size of the grid
        K_original(:,i) = reshape(field1(1:Grid_N(1), 1:Grid_N(2), 1:Grid_N(3)), Grid_N(1)*Grid_N(2)*Grid_N(3), 1);
    else
        % Homogeneous Permeability
        K_original(:,i) = K_original(:,1)*Perm(i);
    end
end