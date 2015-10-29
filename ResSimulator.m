%MATTEO'S RESERVOIR SIMULATOR

%%%%%%%%%%%%%PROPERTIES OF THE RESERVOIR%%%%%%%%%%%%%%%%
Problem = 'Homogeneous';
[Grid, Inj, Prod, Lx, Ly, K]=ReservoirProperties(Problem);
Grid.N = Grid.Nx*Grid.Ny;   %Total number of cells
Grid.Lx = Lx;
Grid.Ly = Ly;

%%%%%%%%%%%%%FLUID PROPERTIES%%%%%%%%%%%%%%%%
fplot = 0;    %if 1 the fractional function curves are plotted
Fluid = FluidProperties(fplot);
%%%%Properties of Injected fluid%%%%
[Inj.Mw, Inj.Mo, Inj.dMw, Inj.dMo] = Mobilities(1,Fluid);

%%%%%%%%%%%%%%%SIMULATOR'S SETTINGS%%%%%%%%%%%
[T, TimeSteps, Options, Tol, Strategy, Sequential, FIM, PlotSolution] =...
    SimulatorSettings();
Grid_Strategy = 'ADM';

%Remove some warnings 
warning('off', 'MATLAB:singularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

%%%%%%%%%%%%%%%INITIAL CONDITIONS%%%%%%%%%%%%%
P = zeros(Grid.Nx, Grid.Ny, 1);
S = ones(Grid.Nx, Grid.Ny, 1)*0.0;
Inj.water = zeros(TimeSteps+1,1);
Prod.water = zeros(TimeSteps+1,1);
Prod.oil = zeros(TimeSteps+1,1);
CumulativeTime = zeros(TimeSteps, 1); 

%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%
t = 0;    %Simulation time
Ndt = 1;  %keeps track of the number of timesteps
Converged = 0;

%%%%%%Choose whether to use ADM or not
switch (Grid_Strategy)
    case ('ADM')
        %ADM Settings
        Pressure_Interpolator =  InterpolationOption; %Constant, Homogeneous, MS
        ADMSettings.maxLevel = 2; %levelMax;
        ADMSettings.tol = 0.2;%tolerance;
        Coarsening = [3, 3; 9, 9; 27, 27]; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
        FIM.ActiveCells = zeros(TimeSteps, ADMSettings.maxLevel + 1);
        %Directory = strcat('../Output/');
        disp('ADM SIMULATION');
        TotalStart = tic;
        ADM;
        TotalTime = toc(TotalStart);
    case ('BaseGrid')
        disp('Base Grid SIMULATION');
        Directory = strcat('../ResultsAndImages/PressureConstrained/', Problem, '/FIM/'); % Output is saved in this directory
        TotalStart = tic;
        BaseGrid;
        TotalTime = toc(TotalStart);
end
disp(char(10));
disp(['The Total Simulation time is ' num2str(TotalTime) ' s']);