%%% Run simulations %%%
Directory_FS = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_FS/';
Directory_ADM = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_ADM/';
Directory_FS_noPc = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_noPc_FS/';
Directory_ADM_noPc = '../Results/4_Papers_Results/AWR_paper/Case4/AWR_Case4_noPc_ADM/';
File = 'Case4.txt';

% run simulations
DARSim2ResSim(Directory, File);
