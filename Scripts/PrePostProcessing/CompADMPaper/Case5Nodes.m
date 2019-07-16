Directory = '../Results/4_Papers_Results/AWR_paper/Case5/Case5_2/ADM/Output/';

Nodes = load(strcat(Directory, 'ADMStats.txt'));
AV_Nodes = mean(Nodes(:, end));
disp(['Case5 av. nodes: ', num2str(AV_Nodes)]);

Stats = load(strcat(Directory, 'SolverStats.txt'));
NL_ADM = sum(Stats(:, 3));
disp(['Case5 adm NL iter: ', num2str(NL_ADM)]);
disp(['Case5 adm time-steps: ', num2str(Stats(end,1))]);

Directory = '../Results/4_Papers_Results/AWR_paper/Case5/Case5_2/FS/Output/';
Stats = load(strcat(Directory, 'SolverStats.txt'));
NL_FS = sum(Stats(:, 3));
disp(['Case5 fs NL iter: ', num2str(NL_FS)]);
disp(['Case5 fs time-steps: ', num2str(Stats(end,1))]);