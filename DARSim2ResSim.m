%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim2ResSim: Reservoir Simulator main file
%
% Requires Matlab 2016b or newer
%
% Run instructions:
% Run from the folder src
% Use addpath(genpath(pwd)); to add all folders to the path
% DARSim2ResSim('Directory', 'FileName', 'PermeabilityDirectory');
% Example of Directory and File
% ImmDirectory = '../Input/ImmHomo/'
% ImmFile = 'ImmHomo';
% PermDir = '../Permeability/'
function ResSimulator = DARSim2ResSim(Directory, File, PermDirectory)
clc;

% Make sure you are in the correct folder
CurrentDir = pwd();
if ~strcmp(CurrentDir(end-10:end), 'DARSim2/src') && ~strcmp(CurrentDir(end-10:end), 'DARSim2\src')
    error('DARSim2 error: you have to be in the folder src to run the code!');
end
clear CurrentDir

addpath(genpath(pwd))

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
fclose('all');

end
