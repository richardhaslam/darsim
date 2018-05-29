%%% Run simulations %%%
%disp('Run 1D')
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_fs/','ImmHomo.txt');
DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_adm_dfdt/','ImmHomo.txt', '../Permeability/');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Small_dt/Homo1D_adm_dfdx/','ImmHomo.txt');

% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Large_dt/Homo1D_fs/','ImmHomo.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Large_dt/Homo1D_adm_dfdt/','ImmHomo.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_1D/Large_dt/Homo1D_adm_dfdx/','ImmHomo.txt');

% Run 2D
% disp('Run 2D')
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_fs/','ImmHomo.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_adm_dfdt/','ImmHomo.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Small_dt/Homo2D_adm_dfdx/','ImmHomo.txt');

% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Large_dt/Homo2D_fs/','ImmHomo.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Large_dt/Homo2D_adm_dtdt/','ImmHomo.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Homogeneous_2D/Large_dt/Homo2D_adm_dfdx/','ImmHomo.txt');

% DARSim2ResSim('../Input/SATINTERPOLATOR/Heterogeneous_noPc/Small_dt/Hetero2D_fs/','ImmHetero.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Heterogeneous_noPc/Small_dt/Hetero2D_adm_dfdt/','ImmHetero.txt');
% DARSim2ResSim('../Input/SATINTERPOLATOR/Heterogeneous_noPc/Small_dt/Hetero2D_adm_dfdx/','ImmHetero.txt');