%%% Run simulations %%%
Directory {1}= '../Results/4_Papers_Results/AWR_paper/Case3/1_AWR_Case3_Small/';
Directory {2}= '../Results/4_Papers_Results/AWR_paper/Case3/2_AWR_Case3_Medium/';
Directory {3}= '../Results/4_Papers_Results/AWR_paper/Case3/3_AWR_Case3_Large/';
% run simulations
File{1} = 'Case3.txt';
File{2} = 'Case3.txt';
File{3} = 'Case3.txt';


parpool(3);

parfor i=1:3
    DARSim2ResSim(Directory{i}, File{i});
end