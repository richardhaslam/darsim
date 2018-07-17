%%% Run simulations %%%
%disp('Run 1D')
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_fs/','ImmHomo.txt', '../Permeability/');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_adm_dfdt/','ImmHomo.txt', '../Permeability/');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_adm_dfdx/','ImmHomo.txt', '../Permeability/');

% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Large_dt/Homo1D_fs/','ImmHomo.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Large_dt/Homo1D_adm_dfdt/','ImmHomo.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Large_dt/Homo1D_adm_dfdx/','ImmHomo.txt');

% Run 2D
% disp('Run 2D')
DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_fs/','ImmHomo.txt', '../Permeability/');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_adm_dfdt/','ImmHomo.txt', '../Permeability/');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_adm_dfdx/','ImmHomo.txt', '../Permeability/');

% DARSim2ResSim('../Input/SATINTERPOLATOR/Barriers_2D/Barriers_2D_fs/'  ,'Barriers_2D.txt', '../Permeability/');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Barriers_2D/Barriers_2D_dfdx/','Barriers_2D.txt', '../Permeability/');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Barriers_2D/Barriers_2D_dfdt/','Barriers_2D.txt', '../Permeability/');

% DARSim2ResSim('../Input/SATINTERPOLATOR/Heterogeneous_noPc/Small_dt/Hetero2D_fs/','ImmHetero.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Heterogeneous_noPc/Small_dt/Hetero2D_adm_dfdt/','ImmHetero.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Heterogeneous_noPc/Small_dt/Hetero2D_adm_dfdx/','ImmHetero.txt');
DARSim2ResSim('D:/SURFdrive/Simulation/DARSim2/Input_Desktop_Dell/ImmHomo/', 'ImmHomo.txt', []);