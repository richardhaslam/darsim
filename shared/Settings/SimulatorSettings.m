function [T, TimeStep, Options, Tol, Strategy, Sequential, FIM, PlotSolution] = SimulatorSettings()
%SIMULATOR SETTINGS
T = 200*24*3600;          %Total time of the simulation [s]
TimeStep = 100000;
Tol=10^(-3);

Strategy = 'FIM'; %Sequential or FIM

%%%%Sequential strategy settings
Sequential.CFL=0.2;  %CFL number
Sequential.MaxExtIter=10;
Sequential.Tol=Tol;
%Implicit Solver: it's used if implicit saturation is required
Sequential.ImpSat=1; %If 1 implicit transport is used
Sequential.ImplicitSolver.fluxfunction = 2; %If 1 it uses the 2nd derivative of dfdS
Sequential.ImplicitSolver.tol=10^(-6);
Sequential.ImplicitSolver.maxNewton=20;

%%%%FIM settings
FIM.MaxIter = 12;
FIM.Tol = 1e-6;
FIM.CFL = 3;

%%%%Plotting options
PlotSolution = 0; % 0 or 1: if 0 no plot, 1 solution is plot during the simulation
Options.Pressure_3D = 0; % 0 or 1, if 1 pressure plot in 3D
Options.problem_1D = 0; % if 1, the plotting for a 1D problem is activated
Options.ContourPlot = 0; % If 1 dynamic contour plot
Options.PlotResiduals = 0; % If 1 Residuals are plotted

%%%Initialise objects for outputting statistics
Sequential.ImplicitSolver.timestep = 0;
Sequential.ImplicitSolver.Newtons = 0;
Sequential.ImplicitSolver.Chops = 0;
FIM.timestep = zeros(TimeStep,1);
FIM.Iter = zeros(TimeStep,1);
FIM.Chops = zeros(TimeStep,1);
end