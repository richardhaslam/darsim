%Choose directories and files
% Folders FS
Directory = '../Results/4_Papers_Results/AWR_paper/Case5/';
for j=1:20
    if ~exist(strcat(Directory,'Case5_',num2str(j)), 'dir');
        mkdir(strcat(Directory,'Case5_',num2str(j)), 'FS');
        mkdir(strcat(Directory,'Case5_',num2str(j)), 'ADM');
    end
    if ~strcmp(strcat('Case5_', num2str(j)), 'Case5_1')
        copyfile(strcat(Directory,'Case5_1/FS/Case5.txt'), strcat(Directory, 'Case5_',num2str(j),'/FS'));
        copyfile(strcat(Directory,'Case5_1/FS/Case5.txt'), strcat(Directory, 'Case5_',num2str(j),'/ADM'));
        copyfile(strcat(Directory,'Case5_1/FS/SimulatorSettings.txt'), strcat(Directory, 'Case5_',num2str(j),'/FS'));
        copyfile(strcat(Directory,'Case5_1/ADM/SimulatorSettings.txt'), strcat(Directory, 'Case5_',num2str(j),'/ADM'));
    end
    %% Modify input files FS
    Dir_fs = strcat(Directory, 'Case5_',num2str(j),'/FS');
    % Read txt into cell A
    A = regexp(fileread(strcat(Dir_fs,'/Case5.txt')), '\n', 'split');
    
    % Change Title
    A{1} = sprintf('%s', '-- DARSim research code: INPUT file');
    % Change Title
    tmp = strfind(A, 'TITLE'); % Search a specific string and find all rows containing matches
    index = find(~cellfun('isempty', tmp));
    A{index + 1} = strcat('Case5_fs_',num2str(j));
    % Change Perm field
    tmp = strfind(A, 'PERMX');
    index = find(~cellfun('isempty', tmp));
    A{index + 1} = strcat('EggModel/EggPerm', num2str(j),'.txt');
    
    % Write cell A into txt
    fid = fopen(strcat(Dir_fs,'/Case5.txt'), 'w');
    fprintf(fid, '%s\n', A{:});
    fclose(fid);
    
    %% Modify input files ADM
    Dir_adm = strcat(Directory, 'Case5_',num2str(j),'/ADM');
    % Read txt into cell A
    A = regexp(fileread(strcat(Dir_adm,'/Case5.txt')), '\n', 'split');
    
    % Change Title
    A{1} = sprintf('%s', '-- DARSim research code: INPUT file');
    % Change Title
    tmp = strfind(A, 'TITLE'); % Search a specific string and find all rows containing matches
    index = find(~cellfun('isempty', tmp));
    A{index + 1} = strcat('Case5_adm_',num2str(j));
    % Change Perm field
    tmp = strfind(A, 'PERMX');
    index = find(~cellfun('isempty', tmp));
    A{index + 1} = strcat('EggModel/EggPerm', num2str(j),'.txt');
    
    % Write cell A into txt
    fid = fopen(strcat(Dir_adm,'/Case5.txt'), 'w');
    fprintf(fid, '%s\n', A{:});
    fclose(fid);
    
end

