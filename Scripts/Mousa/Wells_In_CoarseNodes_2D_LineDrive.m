clc;

%% Number of Injetors/Producers ('QuarterFiveSpot' , 'LineDriveTwo' , 'LineDriveFull')
Well_Pattern = 'QuarterFiveSpot'; % injectors at left and producers at right

Directory = '../Input/SinglePhase/';
File_InputSettings = 'SinglePhase_Original.txt';
File_SimulSettings = 'SimulatorSettings.txt';

File_InputSettings = strcat(Directory,'/',File_InputSettings);
File_SimulSettings = strcat(Directory,'/',File_SimulSettings);

%% Reading SimulatorSettings
fileID = fopen(File_SimulSettings, 'r');
SimulSettingsMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

%% Reading InputSettings
fileID = fopen(File_InputSettings, 'r');
InputSettingsMatrix = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

%% Checking MMs Status
temp = strfind(SimulSettingsMatrix{1}, 'MMs');
MMs = find(~cellfun('isempty', temp));
if str2double(SimulSettingsMatrix{1}(MMs + 1)) == 1
    % Collecting MMs settings
    temp = strfind(SimulSettingsMatrix{1}, 'LEVELS');
    x = find(~cellfun('isempty', temp));
    MMsLevel = str2double(SimulSettingsMatrix{1}(x+1));
    obj.MMsSettings.Coarsening = zeros(3,MMsLevel);
    temp = strfind(SimulSettingsMatrix{1}, 'COARSENING_RATIOS');
    x = find(~cellfun('isempty', temp));
    cx = str2double(SimulSettingsMatrix{1}(x+1))^MMsLevel;
    cy = str2double(SimulSettingsMatrix{1}(x+2))^MMsLevel;
    cz = str2double(SimulSettingsMatrix{1}(x+3))^MMsLevel;
    
    % Collecting Grid size data
    temp = strfind(InputSettingsMatrix{1}, 'SPECGRID');
    x = find(~cellfun('isempty', temp));
    Nx = str2double(InputSettingsMatrix{1}(x+1));
    Ny = str2double(InputSettingsMatrix{1}(x+2));
    Nz = str2double(InputSettingsMatrix{1}(x+3));
    Ncx = Nx/cx;
    Ncy = Ny/cy;
    Ncz = Nz/cz;
    
    % Collecting well properties
    temp = strfind(InputSettingsMatrix{1}, 'INJ1');
    x_inj = find(~cellfun('isempty', temp));
    temp = strfind(InputSettingsMatrix{1}, 'PROD1');
    x_prod = find(~cellfun('isempty', temp));

    % Creating New Input File
    delete(strcat(Directory,'/SinglePhase.txt'));
    fileID = fopen(strcat(Directory,'/SinglePhase.txt'),'w');
    
    % Locating the line where wells will be written
    temp = strfind(InputSettingsMatrix{1}, '--Wells');
    x = find(~cellfun('isempty', temp));
    fprintf(fileID,'%s\n',InputSettingsMatrix{1}{1:x});
    
    switch Well_Pattern
        case('QuarterFiveSpot')
            Well_List_inj  = [ 1   ];
            Well_List_prod = [ Ncy ];
        case('LineDriveTwo')
            Well_List_inj  = [ 1  Ncy ];
            Well_List_prod = [ 1  Ncy ];
        case('LineDriveFull')
            Well_List_inj  = [ 1:Ncy ];
            Well_List_prod = [ 1:Ncy ];
    end
    
    % Writing injection wells
    for count = 1:length(Well_List_inj)
        n = Well_List_inj(count);
        fprintf(fileID,'%s\n',strcat('INJ',num2str(count)));
        i_inj = round(cx/2);
        j_inj = round(cy/2) + (n-1)*cy;
        fprintf(fileID,'%s\n%s\n',num2str(i_inj),num2str(i_inj));
        fprintf(fileID,'%s\n%s\n',num2str(j_inj),num2str(j_inj));
        fprintf(fileID,'%s\n%s\n',num2str(1    ),'NZ');
        fprintf(fileID,'%s\n'    ,InputSettingsMatrix{1}{x_inj+7});
        fprintf(fileID,'%s\n'    ,InputSettingsMatrix{1}{x_inj+8});
        fprintf(fileID,'%s\n'    ,InputSettingsMatrix{1}{x_inj+9});
        fprintf(fileID,'%s\n\n'  ,InputSettingsMatrix{1}{x_inj+10});
    end
    
    % Writing production wells
    for count = 1:length(Well_List_prod)
        n = Well_List_prod(count);
        fprintf(fileID,'%s\n',strcat('PROD',num2str(count)));
        i_prod = Nx - round(cx/2) +1;
        j_prod = round(cy/2) + (n-1)*cy;
        fprintf(fileID,'%s\n%s\n',num2str(i_prod),num2str(i_prod));
        fprintf(fileID,'%s\n%s\n',num2str(j_prod),num2str(j_prod));
        fprintf(fileID,'%s\n%s\n',num2str(1    ),'NZ');
        fprintf(fileID,'%s\n'    ,InputSettingsMatrix{1}{x_prod+7});
        fprintf(fileID,'%s\n'    ,InputSettingsMatrix{1}{x_prod+8});
        fprintf(fileID,'%s\n'    ,InputSettingsMatrix{1}{x_prod+9});
        fprintf(fileID,'%s\n\n'  ,InputSettingsMatrix{1}{x_prod+10});
    end
    
    temp = strfind(InputSettingsMatrix{1}, '--Fractures');
    x = find(~cellfun('isempty', temp));
    fprintf(fileID,'%s\n',InputSettingsMatrix{1}{x:end});
    fclose(fileID);
end
fprintf('Successful!\n');