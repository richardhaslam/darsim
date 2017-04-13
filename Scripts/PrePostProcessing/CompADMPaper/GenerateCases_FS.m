%Choose directories and files
% Folders FS
Directory = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS_natural/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];

for i = 1:5
    for j=1:20
        if ~exist(strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j)), 'dir');
            mkdir(strcat(Directory,Angle(i,:),'/'), strcat(Angle(i,:),'_',num2str(j)));
        end
        if ~strcmp(strcat(Angle(i,:),'_',num2str(j)), '00deg_1')
            copyfile(strcat(Directory,'00deg/00deg_1/BOHetero.txt'), strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j)));
            copyfile(strcat(Directory,'00deg/00deg_1/SimulatorSettings.txt'), strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j)));
        end
        %% Modify input files
        Dir = strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j));
        % Read txt into cell A
        A = regexp( fileread(strcat(Dir,'/BOHetero.txt')), '\n', 'split');
        
        % Change Title
        
        A{1} = sprintf('%s', '-- DARSim research code: INPUT file');
        % Change Title
        tmp = strfind(A, 'TITLE'); % Search a specific string and find all rows containing matches
        index = find(~cellfun('isempty', tmp));
        A{index + 1} = strcat('Case1_',Angle(i,:),'_', num2str(j),'_FS_natural');
        % Change Perm field
        tmp = strfind(A, 'PERMX');
        index = find(~cellfun('isempty', tmp));
        A{index + 1} = strcat(Angle(i,:),'/','Perm_',Angle(i,:),'_', num2str(j), '.txt');
        
        % Write cell A into txt
        fid = fopen(strcat(Dir,'/BOHetero.txt'), 'w');
        fprintf(fid, '%s\n', A{:});
        fclose(fid);
        
        
    end
end


