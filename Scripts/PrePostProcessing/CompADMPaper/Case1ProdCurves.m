%% Compare production curves
close all

Directory_FS = '../Results/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS/';
Directory_ADM_ms = '../Results/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_MS/';
Directory_ADM_const = '../Results/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_Constant/';

PlotProductionCurves(Directory_FS, 'red', 1);
PlotProductionCurves(Directory_ADM_ms, 'blue', 1);
PlotProductionCurves(Directory_ADM_const, 'green', 1);
%% Rescale for latex plot
Dir_fs = strcat(Directory_FS, '/Output/');
Ph1_prod_fs = load(strcat(Dir_fs, 'Prod_Phase1.txt'));
Ph2_prod_fs = load(strcat(Dir_fs, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case1_Prod_Ph1_fs.txt'), 'w');
time = Ph1_prod_fs(:, 1);
index = find(time == 18);
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_fs(1:index, 1), Ph1_prod_fs(1:index,end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case1_Prod_Ph2_fs.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_fs(1:index,1), Ph2_prod_fs(1:index,end)]');
fclose(fid);
% ADM_const
Dir_adm = strcat(Directory_ADM_const,'/Output/');
Ph1_prod_adm = load(strcat(Dir_adm, 'Prod_Phase1.txt'));
Ph2_prod_adm = load(strcat(Dir_adm, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case1_Prod_Ph1_adm_const.txt'), 'w');
time = Ph1_prod_adm(:, 1);
index = find(time == 18);
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_adm(1:index, 1), Ph1_prod_adm(1:index,end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case1_Prod_Ph2_adm_const.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_adm(1:index,1), Ph2_prod_adm(1:index,end)]');
fclose(fid);

% ADM_ms
Dir_ADM = strcat(Directory_ADM_ms, '/Output/');
Ph1_prod = load(strcat(Dir_ADM, 'Prod_Phase1.txt'));
Ph2_prod = load(strcat(Dir_ADM, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case1_Prod_Ph1_adm_ms.txt'), 'w');
time = Ph1_prod(:, 1);
index = find(time == 18);
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod(1:index, 1), Ph1_prod(1:index,end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case1_Prod_Ph2_adm_ms.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod(1:index,1), Ph2_prod(1:index,end)]');
fclose(fid);