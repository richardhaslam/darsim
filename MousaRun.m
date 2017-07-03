% The Run File for Mousa
cd 'D:\Git\DARSim 2\src';

clc;
Directory = '../Input/SinglePhase/';
File = 'SinglePhase.txt';
Directory = '../Input/ImmHomo/';
File = 'ImmHomo.txt';

DARSim2Res = DARSim2ResSim(Directory, File);