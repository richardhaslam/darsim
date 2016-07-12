%DARSim 2 RESERVOIR SIMULATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 9 April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Status = ResSimulator(InputDirectory, InputFile)
clc;
%cd('/media/matteo/LinuxData/PhD/MatteoResSim/src');
addpath(genpath('../src'));

%Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%%%%%%%%%%%%%%%% READ DATA from INPUT file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Errors
Errors = 0;
delete(strcat(InputDirectory, '/Output/RunDiary.txt')); 
diary(strcat(InputDirectory, '/Output/RunDiary.txt'));
disp('******************************************************************');
disp('********************DARSIM 2 RESERVOIR SIMULATOR********************');
disp('******************************************************************');
disp(char(5));
disp('--------------------------------------------');
disp('Reading input file');

[Output, T, Grid, K, Fluid, Inj, Prod, Coupling,...
    Summary, ADMSettings, FlashSettings, TimeSteps] =...
ReadInputFile(InputDirectory, InputFile);

switch (Errors)
    case(1)
        disp('-----------------------------------------------')
        disp('Run failed: the input file contains errors')
    case(0)
        disp('Input file read without any error');
        disp('-----------------------------------------------')
        disp(char(5));
        disp(['******************', num2str(Output.ProblemName),'******************']);
        disp(['Lx ', num2str(Grid.Lx), ' m']);
        disp(['Ly ', num2str(Grid.Ly), ' m']);
        disp(['h  ', num2str(Grid.h), ' m']);
        disp(['Grid: ', num2str(Grid.Nx), ' x ',  num2str(Grid.Ny), ' x 1']);
        disp(char(5));
        disp('---------Simulator Settings-----------');
        disp(['Coupling: ', Coupling.Name]);
        if ~exist(strcat(InputDirectory,'/Output/VTK/'), 'dir')
            mkdir(InputDirectory,'/Output/VTK');
        else
            delete(strcat(InputDirectory,'/Output/VTK/*.vtk'));
        end
        
        %Cances function if capillarity is used
        if (~isempty(Fluid.Pc) && ~strcmp(Fluid.Pc,'JLeverett'))
            Fluid = ComputeCancesFunction(Fluid, reshape(K(1,:,:), Grid.Nx, Grid.Ny), Grid.por);
        end
        
        %%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%% 
        Status = InitializeReservoir(Fluid, Grid, K, FlashSettings);
        Output.Plotter.PlotInitialStatus(Status, Grid, K, Inj, Prod);
   
        %%%%%%%%%%%%%% ADM SETUP %%%%%%%%%%%%%%%%%%
        if (strcmp(Coupling.Name, 'FIM') == 1 && ADMSettings.active == 1)
            disp('ADM Settings');
            disp(['Pressure Interpolator: ', ADMSettings.Pressure_Interpolator]);
            disp(['Number of levels: ', num2str(ADMSettings.maxLevel)]);
            disp(char(5));
            %disp(['Coarsening ratio: ', num2str(ADMSettings)]);
            [Grid, CoarseGrid] = ADMSetup(Grid, K, ADMSettings, Inj, Prod);
        else
            CoarseGrid = 0;
            ADMSettings.maxLevel = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%% TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%
        TotalStart = tic;
        [Status, Summary] = TimeDriver...
            (Grid, K, Fluid, Status, Inj, Prod, T, TimeSteps, Summary, Coupling,...
             Output, CoarseGrid, ADMSettings, FlashSettings);
        TotalTime = toc(TotalStart);
        
        %%%%%%%%%%%%%%%%%%% OUTPUT STATS %%%%%%%%%%%%%%%%%%%%%%
        Output.WriteSummary(Summary);
                    
        %%%%%%%%%%%%%%%%%% END of SIMULATION %%%%%%%%%%%%%%%%%
        disp(char(10));
        disp(['The Total Simulation time is ' num2str(TotalTime) ' s']);
end
diary off
end