% The Run File for Mousa
close all; clc;
fclose('all');
Directory = '../FracGen_IO/';
File = 'Fracture_Input.txt';
Fracture_Network = DARSim2FracGen(Directory, File);


fprintf('******************************************************************\n');

% Destination_Dir = '../Input/SinglePhase/';
% disp(['Copying the output file into ', Destination_Dir]);
% copyfile('../FracGen_IO/Fracture_Output.txt',Destination_Dir);
% disp('The file is copied successfully.');

Destination_Dir = '../Input/ImmHomo/';
disp(['Copying the output file into ', Destination_Dir]);
copyfile('../FracGen_IO/Fracture_Output.txt',Destination_Dir);
disp('The file is copied successfully.');