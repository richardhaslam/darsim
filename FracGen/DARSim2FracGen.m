% Fracture_Generator Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Mousa HosseiniMehr
%TU Delft
%Created: 2017-03-10
%Last modified: 2017-03-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FracGenerator = DARSim2FracGen(Directory, File)
clc;
% Example of Directory and File
% ImmDirectory = '../Input/ImmHomo'
% ImmFile = 'ImmHomo';

% Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%% Set up Diary
if exist(   strcat(Directory, '/Output/Fracture_Output.txt'), 'file')
    delete( strcat(Directory, '/Output/Fracture_Output.txt')        );
end

%% Print title
disp('********************************************************************************');
disp('******************** Fracture Generator for DARSim2 Simulator ******************');
disp('********************************************************************************');
disp(char(5));
disp(['Reading input file ', File, ' from ', Directory]);
disp(char(5));

%% Build objects
% Build Simulator
FracGenerator = Fracture_Generator(Directory, File);
% Read Input File
FracGenerator.Reader.ReadInputFile();
% Build objects
FracGenerator.BuildObjects();
% Print info to screen
FracGenerator.PrintInfo();

%% Run Simulation
TotalStart = tic;
FracGenerator.Run();
TotalTime = toc(TotalStart);

%% Display elapsed time
disp(['The total computation time is ' num2str(TotalTime) ' s']); fprintf('\n');

%% Output Results
FracGenerator.OutputResults();

end