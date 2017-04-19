% DARSim 2 Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 28 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ResSimulator = DARSim2ResSim(Directory, File)
clc;
addpath(genpath(pwd));
% Example of Directory and File
%ImmDirectory = '../Input/ImmHomo'
%ImmFile = 'ImmHomo';

%Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%% Set up Diary
if ~exist(strcat(Directory, '/Output/'), 'dir')
    mkdir(strcat(Directory,'/'), 'Output');
else
    delete(strcat(Directory, '/Output/*.txt'));
end
delete(strcat(Directory, '/Output/RunDiary.txt'));
diary(strcat(Directory, '/Output/RunDiary.txt'));

%% Print title
disp('******************************************************************');
disp('********************DARSIM 2 RESERVOIR SIMULATOR******************');
disp('******************************************************************');
disp(char(5));
disp(['Reading input file ', File, ' from ', Directory]);
disp(char(5));

%% Build objects
% Build Simulator
ResSimulator = Reservoir_Simulator(Directory, File);
% Read Input File
ResSimulator.Reader.ReadInputFile();
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
disp(char(10));
disp(['The Total Simulation time is ' num2str(TotalTime) ' s']);

%Tun off diary
diary off
end