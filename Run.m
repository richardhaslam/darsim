%%% Run simulations %%%
nf = 2;
Directory = cell(nf, 1);
File = cell(nf, 1);
% Choose directories and files
Directory{2} = '../ResultsAndImages/4_Papers_Results/RSC_CompPaper/TestCase2/BOHetero_FS';
File{2} = 'BOHetero.txt';
Directory{1} = '../ResultsAndImages/4_Papers_Results/RSC_CompPaper/TestCase2/BOHetero_ADM';
File{1} = 'BOHetero.txt';

% run simulations
for i = 1:nf
     DARSim2ResSim(Directory{i}, File{i});
 end

%DARSim2ResSim('../Input/BOHomo', 'BOHomo.txt');