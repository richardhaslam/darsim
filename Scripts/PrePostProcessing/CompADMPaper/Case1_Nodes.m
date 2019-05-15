Directory_Small = '../Results/4_Papers_Results/AWR_paper/Case2/1_AWR_Case2_Small/Output/';
Directory_Medium = '../Results/4_Papers_Results/AWR_paper/Case2/2_AWR_Case2_Medium/Output/';
Directory_Large = '../Results/4_Papers_Results/AWR_paper/Case2/3_AWR_Case2_Large/Output/';

%% Rescale for latex plot
Inj_Small = load(strcat(Directory_Small, 'Inj_Phase1.txt'));
Inj_Small = Inj_Small(:, end) ./ (0.5 * Inj_Small(end,end));

Inj_Medium = load(strcat(Directory_Medium, 'Inj_Phase1.txt'));
Inj_Medium = Inj_Medium(:, end) ./ (0.5 * Inj_Medium(end,end));

Inj_Large = load(strcat(Directory_Large, 'Inj_Phase1.txt'));
Inj_Large = Inj_Large(:, end) ./ (0.5 * Inj_Large(end,end));

Nodes_Small = load(strcat(Directory_Small, 'ADMStats_Small.txt'));
Nodes_Small = Nodes_Small(:, end);

Nodes_Medium = load(strcat(Directory_Medium, 'ADMStats_Medium.txt'));
Nodes_Medium = Nodes_Medium(:, end);

Nodes_Large = load(strcat(Directory_Large, 'ADMStats_Large.txt'));
Nodes_Large = Nodes_Large(:, end);

fid = fopen(strcat('Nodes_Small.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Inj_Small, Nodes_Small]');
fclose(fid);

fid = fopen(strcat('Nodes_Medium.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Inj_Medium, Nodes_Medium]');
fclose(fid);

fid = fopen(strcat('Nodes_Large.txt'), 'w');
fprintf(fid, '%10.4f %10.4f\n', [Inj_Large, Nodes_Large]');
fclose(fid);