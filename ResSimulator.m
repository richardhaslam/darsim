%MATTEO'S RESERVOIR SIMULATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%READ DATA from INPUT file%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd('../Code')
InputDirectory = '../Input/Homo_Cap/';
InputFile = strcat(InputDirectory, '/Homogeneous.txt');
ReadInputFile;
if ~exist(strcat(InputDirectory,'/Output/VTK/'), 'dir')
  mkdir(InputDirectory,'/Output/VTK');
end
Directory = strcat(InputDirectory,'/Output/');

%%Plot Permeability Field
%K(1:end/2) = K(1:end/2) *1e1;
PlotPermeability(K, Grid);

%Cances function if capillarity is used
if (~isempty(Fluid.Pc))
    Fluid = ComputeCancesFunction(Fluid);
end
%%%%%%%%%%%%%%%INITIAL CONDITIONS%%%%%%%%%%%%%
P = zeros(Grid.Nx, Grid.Ny, 1);
S = ones(Grid.Nx, Grid.Ny, 1)*0.1;

CumulativeTime = zeros(TimeSteps, 1); 

%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%
t = 0;    %Simulation time
Ndt = 1;  %keeps track of the number of timesteps
Converged = 0;
%Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%%%%%%Choose whether to use ADM or not
switch (ADMSettings.active)
    case (1)
        %ADM Settings
        FIM.ActiveCells = zeros(TimeSteps, ADMSettings.maxLevel + 1);
        disp('ADM SIMULATION');
        TotalStart = tic;
        ADM;
        TotalTime = toc(TotalStart);
    case (0)
        disp('Base Grid SIMULATION');
        TotalStart = tic;
        BaseGrid;
        TotalTime = toc(TotalStart);
end
disp(char(10));
disp(['The Total Simulation time is ' num2str(TotalTime) ' s']);