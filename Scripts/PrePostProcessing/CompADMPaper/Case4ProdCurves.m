% Compare production curves
close all

Directory_FS = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_FS/';
Directory_ADM = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_ADM/';
Directory_ADM_dz2 = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_ADM_dz2/';

Directory_noPc_FS = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_noPc_FS/';
Directory_noPc_ADM = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_noPc_ADM/';
Directory_noPc_ADM_dz2 = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_noPc_ADM_dz2/';

% PlotProductionCurves(Directory_FS, 'red');
% PlotProductionCurves(Directory_ADM, 'yellow');
PlotProductionCurves(Directory_noPc_FS, 'blue');
PlotProductionCurves(Directory_noPc_ADM, 'green');
PlotProductionCurves(Directory_noPc_ADM_dz2, 'red');
%% Rescale for latex plot
Dir_fs = strcat(Directory_FS, '/Output/');
Ph1_prod_fs = load(strcat(Dir_fs, 'Prod_Phase1.txt'));
Ph2_prod_fs = load(strcat(Dir_fs, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case4_Prod_Ph1_fs.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_fs(:, 1), Ph1_prod_fs(:,end)/Ph1_prod_fs(end, end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case4_Prod_Ph2_fs.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_fs(:,1), Ph2_prod_fs(:,end)/Ph2_prod_fs(end, end)]');
fclose(fid);
% ADM
Dir_adm = strcat(Directory_ADM,'/Output/');
Ph1_prod_adm = load(strcat(Dir_adm, 'Prod_Phase1.txt'));
Ph2_prod_adm = load(strcat(Dir_adm, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case4_Prod_Ph1_adm_const.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_adm(:, 1), Ph1_prod_adm(:,end)/Ph1_prod_fs(end, end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case4_Prod_Ph2_adm_const.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_adm(:,1), Ph2_prod_adm(:,end/Ph2_prod_fs(end, end))]');
fclose(fid);

% ADM_ms
Dir_ADM = strcat(Directory_ADM_dz2, '/Output/');
Ph1_prod_dz2 = load(strcat(Dir_ADM, 'Prod_Phase1.txt'));
Ph2_prod_dz2 = load(strcat(Dir_ADM, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case4_Prod_Ph1_adm_ms.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_dz2(:, 1), Ph1_prod_dz2(:,end)/Ph1_prod_fs(end, end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case4_Prod_Ph2_adm_ms.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_dz2(:,1), Ph2_prod_dz2(:,end)/Ph2_prod_fs(end, end)]');
fclose(fid);

%% Rescale for latex plot
Dir_fs = strcat(Directory_noPc_FS, '/Output/');
Ph1_prod_fs = load(strcat(Dir_fs, 'Prod_Phase1.txt'));
Ph2_prod_fs = load(strcat(Dir_fs, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case4_noPc_Prod_Ph1_fs.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_fs(:, 1), Ph1_prod_fs(:,end)/Ph1_prod_fs(end, end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case4_noPc_Prod_Ph2_fs.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_fs(:,1), Ph2_prod_fs(:,end)/Ph2_prod_fs(end, end)]');
fclose(fid);
% ADM
Dir_adm = strcat(Directory_noPc_ADM,'/Output/');
Ph1_prod_adm = load(strcat(Dir_adm, 'Prod_Phase1.txt'));
Ph2_prod_adm = load(strcat(Dir_adm, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case4_noPc_Prod_Ph1_adm_const.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_adm(:, 1), Ph1_prod_adm(:,end)/Ph1_prod_fs(end, end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case4_noPc_Prod_Ph2_adm_const.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_adm(:,1), Ph2_prod_adm(:,end/Ph2_prod_fs(end, end))]');
fclose(fid);

% ADM_ms
Dir_ADM = strcat(Directory_noPc_ADM_dz2, '/Output/');
Ph1_prod_dz2 = load(strcat(Dir_ADM, 'Prod_Phase1.txt'));
Ph2_prod_dz2 = load(strcat(Dir_ADM, 'Prod_Phase2.txt'));
% Phase 1
fid = fopen(strcat('Case4_noPc_Prod_Ph1_adm_ms.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph1_prod_dz2(:, 1), Ph1_prod_dz2(:,end)/Ph1_prod_fs(end, end)]');
fclose(fid);
% Phase 2
fid = fopen(strcat('Case4_noPc_Prod_Ph2_adm_ms.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Ph2_prod_dz2(:,1), Ph2_prod_dz2(:,end)/Ph2_prod_fs(end, end)]');
fclose(fid);