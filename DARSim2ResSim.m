% DARSim 2 Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 28 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ResSimulator = DARSim2ResSim(Directory, File, PermDirectory)
clc;
% addpath(genpath(pwd));
% Example of Directory and File
%ImmDirectory = '../Input/ImmHomo'
%ImmFile = 'ImmHomo';

%Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
% close all files if there's any open one
fclose('all');

%% Set up Diary
if ~exist(strcat(Directory, '/Output/'), 'dir')
    mkdir(strcat(Directory,'/'), 'Output');
else
    delete(strcat(Directory, '/Output/*.txt'));
end
diary(strcat(Directory, '/Output/RunDiary.txt'));

%% Print title
disp('******************************************************************');
disp('********************DARSIM 2 RESERVOIR SIMULATOR******************');
disp('******************************************************************');
disp(newline);
disp(['Reading input file ', File, ' from ', Directory]);
disp(newline);

%% Build objects
% Build Simulator
ResSimulator = Reservoir_Simulator(Directory, File, PermDirectory);
% Read Input File
ResSimulator.Reader.ReadInputFile(ResSimulator.Builder);
% Build objects
ResSimulator.BuildObjects();
% Print info to screen
ResSimulator.PrintInfo();

%% Run Simulation
TotalStart = tic;
ResSimulator.Run();
TotalTime = toc(TotalStart);

%% Output Results
ResSimulator.OutputResults();

%% Display elapsed time
disp(newline);
disp(['The Total Simulation time is ' num2str(TotalTime) ' s']);

%Tun off diary
diary off
end