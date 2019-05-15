% Number of Grid cells
Directory_ADM = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_molar/';
Directory_FS = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];
dz = ['0.03';'0.07';'0.15'];
Error_norm = zeros(5, 2, 3);
for i = 1:1
    for j=1:1
        for k=1:3
            %% Directory ADM Solutions
            Dir_adm = strcat(Directory_ADM,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k),'/Output/');
            Name_adm = 'ADMStats.txt';
            Stats = load(strcat(Dir_adm,Name_adm));
            N_blocks = Stats(:,4);
            TimeSteps = 1:1:length(N_blocks);
            fid = fopen(strcat('N_cells_',Angle(i,:),'_dz_', num2str(k),'.txt'), 'w');
            fprintf(fid, '%d %3.2f\n', [TimeSteps', N_blocks]');
            fclose(fid);
        end
    end
end