% DARSim 2 Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 14 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ResSimulator = DARSim2ResSim(Directory, File)
%% Set up Diary
delete(strcat(Directory, '/Output/RunDiary.txt')); 
diary(strcat(Directory, '/Output/RunDiary.txt'));

%% Print title
disp('******************************************************************');
disp('********************DARSIM 2 RESERVOIR SIMULATOR********************');
disp('******************************************************************');
disp(char(5));
disp('--------------------------------------------');
disp('Reading input file');

%% Build objects
% Build Simulator
ResSimulator = Reservoir_Simulator(Directory, File);
% Read Input File
ResSimulator.Reader.ReadInputFile();
% Build objects
ResSimulator.BuildObjects();

%% Run Simulation
TotalStart = tic;
ResSimulator.Run();
TotalTime = toc(TotalStart);

%% Output Results
ResSimulator.OuputResults();

%% Display elapsed time
disp(char(10));
disp(['The Total Simulation time is ' num2str(TotalTime) ' s']);

%Tun off diary
diary off
end