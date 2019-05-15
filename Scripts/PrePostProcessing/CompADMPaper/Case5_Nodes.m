Directory = '../Results/4_Papers_Results/AWR_paper/Case5/Case5_2/ADM/Output/';

%% Rescale for latex plot
Time = load(strcat(Directory, 'Inj_Phase1.txt'));
Time = Time(:, 1);

Nodes = load(strcat(Directory, 'ADMStats.txt'));
Nodes = Nodes(:, end);

fid = fopen(strcat('Nodes_Case5.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Time, Nodes]');
fclose(fid);
