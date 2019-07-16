% Production
Directory_ADM = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_molar/';
Directory_FS = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS/';
FinalTimes = [40, 60, 50, 25, 40];
for i = 1:1
    LastTime = FinalTimes(i)*365;
    for j=1:1
        Dir_fs = strcat(Directory_FS,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/Output/');
        Ph1_prod_fs = load(strcat(Dir_fs, 'Prod_Phase1.txt'));
        Ph2_prod_fs = load(strcat(Dir_fs, 'Prod_Phase2.txt'));
        LastTimeIndex_fs = find(Ph1_prod_fs(:,1) == LastTime);
        fid = fopen(strcat('Scaled_Prod_Ph1_',Angle(i,:),'_fs.txt'), 'w');
        fprintf(fid, '%10.2f %10.2f\n', [Ph1_prod_fs(1:LastTimeIndex_fs,1)/365, Ph1_prod_fs(1:LastTimeIndex_fs,end)/Ph1_prod_fs(LastTimeIndex_fs, end)]');
        fclose(fid);
        fid = fopen(strcat('Scaled_Prod_Ph2_',Angle(i,:),'_fs.txt'), 'w');
        fprintf(fid, '%10.2f %10.2f\n', [Ph2_prod_fs(1:LastTimeIndex_fs,1)/365, Ph2_prod_fs(1:LastTimeIndex_fs,end)/Ph2_prod_fs(LastTimeIndex_fs, end)]');
        fclose(fid);
        for k=1:3
            %% Directory ADM Solutions
            Dir_adm = strcat(Directory_ADM,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k),'/Output/');    
            Ph1_prod_adm = load(strcat(Dir_adm, 'Prod_Phase1.txt'));
            Ph2_prod_adm = load(strcat(Dir_adm, 'Prod_Phase2.txt'));
            LastTimeIndex = find(Ph1_prod_adm(:,1) == LastTime);
            
            fid = fopen(strcat('Scaled_Prod_Ph1_',Angle(i,:),'_dz_', num2str(k),'.txt'), 'w');
            fprintf(fid, '%10.2f %10.2f\n', [Ph1_prod_adm(1:LastTimeIndex, 1)/365, Ph1_prod_adm(1:LastTimeIndex,end)/Ph1_prod_fs(LastTimeIndex_fs, end)]');
            fclose(fid);
            fid = fopen(strcat('Scaled_Prod_Ph2_',Angle(i,:),'_dz_', num2str(k),'.txt'), 'w');
            fprintf(fid, '%10.2f %10.2f\n', [Ph2_prod_adm(1:LastTimeIndex, 1)/365, Ph2_prod_adm(1:LastTimeIndex,end)/Ph2_prod_fs(LastTimeIndex_fs, end)]');
            fclose(fid);
        end
    end
end