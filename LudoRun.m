clc;


% %Homogeneous 2D permeability fiead
% Directory = '../Input/Examples/ImmiscibleTwoPhaseFlow/Homogeneous2D';
% File = 'Homo2D.txt';
% % %  
% DARSim2Res = DARSim2ResSim(Directory, File,'');
% 
% % %Homogeneous 3D permeability fiead
% Directory = '../Input/Examples/ImmiscibleTwoPhaseFlow/Homogeneous3D';
% File = 'Homo3D.txt';
%  
% DARSim2Res = DARSim2ResSim(Directory, File,'');

% % %Homogeneous 2D with barriers permeability fiead
% Directory = '../Input/Examples/ImmiscibleTwoPhaseFlow/Barrier2D';
% File = 'Barrier.txt';

% DARSim2Res = DARSim2ResSim(Directory, File,'../Permeability');


% %SPE10 Top permeability fiead
Directory = '../Input/Examples/ImmiscibleTwoPhaseFlow/SPE10TOP/ADM_fine';
File = 'SPE10TOP.txt';
  
DARSim2Res = DARSim2ResSim(Directory, File,'../Permeability');
