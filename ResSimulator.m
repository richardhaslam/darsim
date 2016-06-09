%MATTEO'S RESERVOIR SIMULATOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 9 April 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ResSimulator(InputDirectory, InputFile) %#ok<INUSD>
clc;
%cd('/media/matteo/LinuxData/PhD/MatteoResSim/src');
addpath(genpath('../src'));

%Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%%%%%%%%%%%%%%%% READ DATA from INPUT file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Errors
Errors = 0;
%InputDirectory = '../Input/GeoStat_ADM';
%InputFile = strcat(InputDirectory, '/GeoStat.txt');
delete(strcat(InputDirectory, '/Output/RunDiary.txt')); 
diary(strcat(InputDirectory, '/Output/RunDiary.txt'));
disp('******************************************************************');
disp('********************DARSIM 2 RESERVOIR SIMULATOR********************');
disp('******************************************************************');
disp(char(5));
disp('--------------------------------------------')
disp('Reading input file')
ReadInputFile;
switch (Errors)
    case(1)
        disp('-----------------------------------------------')
        disp('Run failed: the input file contains errors')
    case(0)
        disp('Input file read without any error');
        disp('-----------------------------------------------')
        disp(char(5));
        disp(['******************', num2str(Problem),'******************']);
        disp(['Lx ', num2str(Grid.Lx), ' m']);
        disp(['Ly ', num2str(Grid.Ly), ' m']);
        disp(['h  ', num2str(Grid.h), ' m']);
        disp(['Grid: ', num2str(Grid.Nx), ' x ',  num2str(Grid.Ny), ' x 1']);
        disp(char(5));
        if ~exist(strcat(InputDirectory,'/Output/VTK/'), 'dir')
            mkdir(InputDirectory,'/Output/VTK');
        end
        Directory = strcat(InputDirectory,'/Output/');
        
        %%Plot Permeability Field
        PlotPermeability(K, Grid);
        
        %Cances function if capillarity is used
        if (~isempty(Fluid.Pc) && ~strcmp(Fluid.Pc,'JLeverett'))
            Fluid = ComputeCancesFunction(Fluid, reshape(K(1,:,:), Grid.Nx, Grid.Ny), Grid.por);
        end
        
        %%%%%%%%%%%%%%% INITIAL CONDITIONS %%%%%%%%%%%%%
        P = zeros(Grid.Nx, Grid.Ny, 1);
        S = ones(Grid.Nx, Grid.Ny, 1)*0.1;
        
        [Status] = Initialize(Fluid,Grid,FlashSettings);
        P = reshape(Status.p,Grid.Nx,Grid.Ny);         %Finds pressure in rectangular matrix form used for plotting and convection
        S = reshape(Status.s(:,1),Grid.Nx,Grid.Ny);    %Finds saturation in rectangular matrix form used for plotting and convection
   
        %%%%%%%%%%%%%% ADM SETUP %%%%%%%%%%%%%%%%%%
        if (strcmp(Strategy, 'FIM') == 1 && ADMSettings.active == 1)
            ADMSetup;
        else
            CoarseGrid = 0;
            maxLevel = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%% TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%
        TotalStart = tic;
        TimeDriver;
        TotalTime = toc(TotalStart);
        
        %%%%%%%%%%%%%%%%%%% OUTPUT STATS %%%%%%%%%%%%%%%%%%%%%%
        OutputStatistics;
        
        %%%%%%%%%%%%%%%%%% END of SIMULATION %%%%%%%%%%%%%%%%%
        disp(char(10));
        disp(['The Total Simulation time is ' num2str(TotalTime) ' s']);
end
diary off
end