% The Run File for Mousa
close all; clc;
fclose('all');
Directory = '../FracGen_IO/';
File = 'Fracture_Input.txt';
Fracture_Network = DARSim2FracGen(Directory, File);


fprintf('******************************************************************\n');
%Destination_Dir = 'D:\Git\DARSim 2\Input\SinglePhase\';
%disp(['Copying the output file into ', Destination_Dir]);
%copyfile 'Output/Fracture_Output.txt' 'D:\Git\DARSim 2\Input\SinglePhase\';
copyfile '../FracGen_IO/Fracture_Output.txt' '../Input/ImmHomo/';
disp(['The file is copied successfully.']);