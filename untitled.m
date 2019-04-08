% DARSim2ResSim('../Input/SATINTERPOLATOR/SPE10B/SPE10B_fs', 'SPE10B.txt', '../Permeability/');
 %DARSim2ResSim('../Input/SATINTERPOLATOR/SPE10B/SPE10B_dfdx', 'SPE10B.txt', '../Permeability/');
% DARSim2FracGen('../FracGen_IO', 'Fracture_Input.txt');
% 
% !cp ../FracGen_IO/Fracture_Output.txt ../Input/ImmHomo

DARSim2ResSim('../Input/ImmHomo', 'ImmHomo.txt', '../Permeability/');