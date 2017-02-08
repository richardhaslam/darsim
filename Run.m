%%% Run simulations %%%
Directory = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_ADM_molar/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];
dz = ['0.07';'0.15';'0.20'];
% run simulations
for i = 1:1
    for j=1:20
        for k=1:3
        Dir = strcat(Directory,Angle(i,:),'/',Angle(i,:),'_',num2str(j),'/dz_',num2str(k));
        File = 'BOHetero.txt';
        DARSim2ResSim(Dir, File);
        end
    end
 end
