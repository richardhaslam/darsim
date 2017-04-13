%% Compare production curves
close all

Directory_FS = '../Results/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS_molar/';
Directory_ADM = '../Results/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_molar/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];

for i = 3:3
    for j=1:20
        Case1 = strcat(Directory_FS, Angle(i,:),'/',Angle(i,:),'_',num2str(j));
        PlotProductionCurves(Case1, 'red');
        for k=1:1
            Case2 = strcat(Directory_ADM, Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k));
            PlotProductionCurves(Case2, 'blue');
        end
    end
end
