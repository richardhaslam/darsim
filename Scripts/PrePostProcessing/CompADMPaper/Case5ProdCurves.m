% Compare production curves
close all

Directory_FS = '../Results/4_Papers_Results/AWR_paper/Case5/Case5_2/FS';
Directory_ADM = '../Results/4_Papers_Results/AWR_paper/Case5/Case5_2/ADM';

PlotProductionCurves(Directory_FS, 'red:', 1);
PlotProductionCurves(Directory_ADM, 'black:', 1);

%% Rescale for latex plot
Dir_fs = strcat(Directory_FS, '/Output/');
Ph1_prod_fs = load(strcat(Dir_fs, 'Prod_Phase1.txt'));
Ph2_prod_fs = load(strcat(Dir_fs, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case5_Prod_Ph1_fs.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_fs(:, 1), Ph1_prod_fs(:,end)/Ph1_prod_fs(end, end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case5_Prod_Ph2_fs.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_fs(:,1), Ph2_prod_fs(:,end)/Ph2_prod_fs(end, end)]');
fclose(fid);
% ADM
Dir_adm = strcat(Directory_ADM,'/Output/');
Ph1_prod_adm = load(strcat(Dir_adm, 'Prod_Phase1.txt'));
Ph2_prod_adm = load(strcat(Dir_adm, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case5_Prod_Ph1_adm.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_adm(:, 1), Ph1_prod_adm(:,end)/Ph1_prod_fs(end, end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case5_Prod_Ph2_adm.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_adm(:,1), Ph2_prod_adm(:,end)/Ph2_prod_fs(end, end)]');
