% Folders Molar ADM
Directory = '../../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_natural/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];
dz = ['0.03';'0.07';'0.15'];
for i = 1:5
    for j=1:20
        for k=1:3
            if ~exist(strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k)), 'dir');
                mkdir(strcat(Directory,Angle(i,:),'/'), strcat(Angle(i,:),'_',num2str(j), '/dz_', num2str(k)));
            end
            if ~strcmp(strcat(Angle(i,:),'_',num2str(j),'_dz_', num2str(k)), '00deg_1_dz_1')
                copyfile(strcat(Directory,'00deg/00deg_1/dz_1/BOHetero.txt'), strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k)));
                copyfile(strcat(Directory,'00deg/00deg_1/dz_1/SimulatorSettings.txt'), strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k)));
            end
            %% Modify input files
            Dir = strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k));
            % Read txt into cell A
            A = regexp( fileread(strcat(Dir,'/BOHetero.txt')), '\n', 'split');
            B = regexp( fileread(strcat(Dir,'/SimulatorSettings.txt')), '\n', 'split');
            
            % Change Title
            
            A{1} = sprintf('%s', '-- DARSim research code: INPUT file');
            % Change Title
            tmp = strfind(A, 'TITLE'); % Search a specific string and find all rows containing matches
            index = find(~cellfun('isempty', tmp));
            A{index + 1} = strcat('Case1_',Angle(i,:),'_', num2str(j),'natural_ADM_dz_', num2str(k));
            % Change Perm field
            tmp = strfind(A, 'PERMX');
            index = find(~cellfun('isempty', tmp));
            A{index + 1} = strcat(Angle(i,:),'/','Perm_',Angle(i,:),'_', num2str(j), '.txt');
            
            % 
            tmp = strfind(B, 'ADM');
            index = find(~cellfun('isempty', tmp));
            B{index + 1} = '1';
            
            tmp = strfind(B, 'TOLERANCE');
            index = find(~cellfun('isempty', tmp));
            B{index + 1} = dz(k,:);
            
            % Write cell A into txt
            fid = fopen(strcat(Dir,'/BOHetero.txt'), 'w');
            fprintf(fid, '%s\n', A{:});
            fclose(fid);
            % Write cell B into txt
            fid = fopen(strcat(Dir,'/SimulatorSettings.txt'), 'w');
            fprintf(fid, '%s\n', B{:});
            fclose(fid);
        end
    end
end