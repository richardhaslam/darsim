%MATTEO'S RESERVOIR SIMULATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 21 March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%%%%%%%%%%%%%%%% READ DATA from INPUT file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd('../Code')
InputDirectory = '../Input/Homogeneous/';
InputFile = strcat(InputDirectory, '/Homogeneous.txt');
ReadInputFile;
if ~exist(strcat(InputDirectory,'/Output/VTK/'), 'dir')
  mkdir(InputDirectory,'/Output/VTK');
end
Directory = strcat(InputDirectory,'/Output/');

%%Plot Permeability Field
PlotPermeability(K, Grid);

%Cances function if capillarity is used
if (~isempty(Fluid.Pc))
    Fluid = ComputeCancesFunction(Fluid, reshape(K(1,:,:), Grid.Nx, Grid.Ny), Grid.por);
end

%%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%
P = zeros(Grid.Nx, Grid.Ny, 1);
S = ones(Grid.Nx, Grid.Ny, 1)*0.1;
CumulativeTime = zeros(TimeSteps, 1);

%%%%%%%%%%%%%% ADM SETUP %%%%%%%%%%%%%%%%%%
if (strcmp(Strategy, 'FIM') == 1 && ADMSettings.active == 1)
    ADMSetup;
else
    CoarseGrid = 0;
end

%%%%%%%%%%%%%%%%%%%%%% TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%
TotalStart = tic;
t = 0;    %Simulation time
Ndt = 1;  %keeps track of the number of timesteps
Converged = 0;
TimeDriver;
TotalTime = toc(TotalStart);

%%%%%%%%%%%%%%%%%%% OUTPUT STATS %%%%%%%%%%%%%%%%%%%%%%
OutputStatistics;

%%%%%%%%%%%%%%%%%% END of SIMULATION %%%%%%%%%%%%%%%%%
disp(char(10));
disp(['The Total Simulation time is ' num2str(TotalTime) ' s']);
