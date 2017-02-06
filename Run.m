%%% Run simulations %%%
Directory = '../ResultsAndImages/4_Papers_Results/AWR_paper/Case1/AWR_Case1_FS/';
Angle = ['00deg'; '15deg'; '30deg'; '45deg'; '98deg'];

% run simulations
for j = 1:20
    for i=1:5
        Dir = strcat(Directory, Angle(i,:),'/',Angle(i,:),'_',num2str(j));
        File = 'BOHetero.txt';
        DARSim2ResSim(Dir, File);
    end
 end
