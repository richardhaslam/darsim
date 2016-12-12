%%% Run simulations %%%
% nf = 4;
% Directory = cell(nf, 1);
% File = cell(nf, 1);
% Choose directories and files
% Directory{3} = '../ResultsAndImages/4_Papers_Results/RSC_CompPaper/TestCase3/BOHetero_FS';
% File{3} = 'BOHetero.txt';
% Directory{4} = '../ResultsAndImages/4_Papers_Results/RSC_CompPaper/TestCase3/BOHetero_ADM';
% File{4} = 'BOHetero.txt';
% 
% Directory{1} = '../ResultsAndImages/4_Papers_Results/RSC_CompPaper/TestCase2/BOHomo_FS';
% File{1} = 'BOHomo.txt';
% Directory{2} = '../ResultsAndImages/4_Papers_Results/RSC_CompPaper/TestCase2/BOHomo_ADM';
% File{2} = 'BOHomo.txt';
% 
% run simulations
% for i = 1:nf
%      DARSim2ResSim(Directory{i}, File{i});
%  end

DARSim2ResSim('../Input/SPE10T', 'SPE10T.txt');