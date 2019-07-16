%%% Run simulations %%%
Directory {1}= '../Results/4_Papers_Results/AWR_paper/Case2/1_AWR_Case2_Small/';
Directory {2}= '../Results/4_Papers_Results/AWR_paper/Case2/2_AWR_Case2_Medium/';
Directory {3}= '../Results/4_Papers_Results/AWR_paper/Case2/3_AWR_Case2_Large/';
% run simulations
File{1} = 'Case2.txt';
File{2} = 'Case2.txt';
File{3} = 'Case2.txt';


parpool(3);

parfor i=1:3
    DARSim2ResSim(Directory{i}, File{i});
end