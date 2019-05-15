Directory_noPc_dz1 = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_noPc_ADM/Output/';
Directory_noPc_dz2 = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_noPc_ADM_dz3/Output/';
Directory_Pc_dz1 = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_ADM/Output/';
Directory_Pc_dz2 = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_ADM_dz3/Output/';

%% Rescale for latex plot
Time_noPc_dz1 = load(strcat(Directory_noPc_dz1, 'Inj_Phase1.txt'));
Time_noPc_dz1 = Time_noPc_dz1(:, 1);

Time_noPc_dz2 = load(strcat(Directory_noPc_dz2, 'Inj_Phase1.txt'));
Time_noPc_dz2 = Time_noPc_dz2(:, 1);

Time_Pc_dz1 = load(strcat(Directory_Pc_dz1, 'Inj_Phase1.txt'));
Time_Pc_dz1 = Time_Pc_dz1(:, 1);

Time_Pc_dz2 = load(strcat(Directory_Pc_dz2, 'Inj_Phase1.txt'));
Time_Pc_dz2 = Time_Pc_dz2(:, 1);

Nodes_noPc_dz1 = load(strcat(Directory_noPc_dz1, 'ADMStats.txt'));
Nodes_noPc_dz1 = Nodes_noPc_dz1(:, end);

Nodes_noPc_dz2 = load(strcat(Directory_noPc_dz2, 'ADMStats.txt'));
Nodes_noPc_dz2 = Nodes_noPc_dz2(:, end);

Nodes_Pc_dz1 = load(strcat(Directory_Pc_dz1, 'ADMStats.txt'));
Nodes_Pc_dz1 = Nodes_Pc_dz1(:, end);

Nodes_Pc_dz2 = load(strcat(Directory_Pc_dz2, 'ADMStats.txt'));
Nodes_Pc_dz2 = Nodes_Pc_dz2(:, end);

fid = fopen(strcat('Nodes_Case4_noPc_dz1.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Time_noPc_dz1, Nodes_noPc_dz1]');
fclose(fid);

fid = fopen(strcat('Nodes_Case4_noPc_dz3.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Time_noPc_dz2, Nodes_noPc_dz2]');
fclose(fid);

fid = fopen(strcat('Nodes_Case4_Pc_dz1.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Time_Pc_dz1, Nodes_Pc_dz1]');
fclose(fid);

fid = fopen(strcat('Nodes_Case4_Pc_dz3.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Time_Pc_dz2, Nodes_Pc_dz2]');
fclose(fid);