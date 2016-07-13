% DARSim 2 Reservoir Simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ResSimulator = DARSim2ResSim(Directory, File)

%% Build objects
% Build Simulator
ResSimulator = Reservoir_Simulator(Directory, File);
% Read Input File
ResSimulator.Reader.ReadInputFile();
%
ResSimulator.BuildObjects();

%% Run Simulation
ResSimulator.Run();

%% Output Results
ResSimulator.OuputResults();

end