Directory = '../Input/ImmHomo_ADM/Output/';

%% Rescale for latex plot
Time = load(strcat(Directory, 'Inj_Phase1.txt'));
Time = Time(:, 1);

Nodes = load(strcat(Directory, 'ADMStats.txt'));
Nodes = Nodes(:, end);

fid = fopen(strcat('Nodes_ImmFractured.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Time, Nodes]');
fclose(fid);