clc;

% %Homogeneous permeability fiead
% Directory = '../ImmHomo/';
% File = 'ImmHomo.txt';
%  
% DARSim2Res = DARSim2ResSim(Directory, File,'');
% 
% Directory = '../Input/SPE10T1/';
% File = 'SPE10T.txt';
%  
% DARSim2Res = DARSim2ResSim(Directory, File,'../Permeability/');

Directory = '../Input/Barrier/';
File = 'Barrier.txt';
 
DARSim2Res = DARSim2ResSim(Directory, File,'../Permeability/');