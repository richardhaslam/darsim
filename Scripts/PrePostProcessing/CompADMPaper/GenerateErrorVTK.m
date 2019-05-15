Directory_FS = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS/';
Directory_ADM = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_molar/';
Directory_comp = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/Comparison/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];
n_files = [40, 60, 50, 25, 40];
for i=1:1
    for j=1:1
        Case1 = strcat(Directory_FS, Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/');
        for k=1:3
            Case2 = strcat(Directory_ADM, Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k),'/');
            ComparisonDir = strcat(Directory_comp,Angle(i,:),'/', Angle(i,:),'_',num2str(j),'/dz_',num2str(k),'/');
            if ~exist(ComparisonDir, 'dir');
                mkdir(ComparisonDir);
            end
            VTKErrorFile(Case1, Case2, ComparisonDir, n_files(i));
        end
    end
end