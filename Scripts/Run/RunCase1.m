%%% Run simulations %%%
Directory_FS = '../Results/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS/';
Directory_ADM_constant = '../Results/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_Constant/';
Directory_ADM_ms = '../Results/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_MS/';

% run simulations
File = 'Case1.txt';
DARSim2ResSim(Directory_FS, File);
