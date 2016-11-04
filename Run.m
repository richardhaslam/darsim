%%% Run simulations %%%
nf = 1;
Directory = cell(nf, 1);
File = cell(nf, 1);
% Choose directories and files
Directory{1} = '../Input/CompHetero';
File{1} = 'CompHetero.txt';
Directory{2} = '../Input/BOHomo';
File{2} = 'BOHomo.txt';

% run simulations
for i = 1:nf
    DARSim2ResSim(Directory{i}, File{i});
end