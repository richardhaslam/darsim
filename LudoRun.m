clc;


% %Homogeneous permeability fiead
% Directory = '../Input/Examples/ImmiscibleTwoPhaseFlow/Homogeneous2D';
% File = 'Homo2D.txt';
% % %  
% DARSim2Res = DARSim2ResSim(Directory, File,'');

% %Homogeneous permeability fiead
Directory = '../Input/Examples/ImmiscibleTwoPhaseFlow/Homogeneous3D';
File = 'Homo3D.txt';
 
DARSim2Res = DARSim2ResSim(Directory, File,'');

% % %Homogeneous permeability fiead
% Directory = '../Input/Examples/ImmiscibleTwoPhaseFlow/Barrier2D';
% File = 'Barrier.txt';
% % %  
% DARSim2Res = DARSim2ResSim(Directory, File,'../Permeability');