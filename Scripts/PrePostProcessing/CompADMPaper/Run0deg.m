%%% Run simulations %%%
Directory_FS = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS/';
Directory_ADM = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_natural/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];
dz = ['0.05';'0.07';'0.10'];
% run simulations
for i = 1:1
    for j=6:6
        parfor k=2:2
            Dir = strcat(Directory_ADM, Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k));
            File = 'BOHetero.txt';
            DARSim2ResSim(Dir, File);
        end
    end
 end
