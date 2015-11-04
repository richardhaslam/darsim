%SIMULATOR SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FIM, Sequential, ADMSettings] = ...
    SimulatorSettings(TimeSteps, Strategy, settings, impsat, adm, inputMatrix)
%SIMULATOR SETTINGS
switch (Strategy)
    case ('Sequential')
        %%%%Sequential strategy settings
        Sequential.MaxExtIter = str2double(inputMatrix{1}(settings + 1));
        Sequential.Tol = str2double(inputMatrix{1}(settings +2));
        Sequential.CFL = str2double(inputMatrix{1}(settings +3));  %CFL number
        if (impsat~=0)
            %Implicit Solver: it's used if implicit saturation is required
            Sequential.ImpSat=1; %If 1 implicit transport is used
            Sequential.ImplicitSolver.fluxfunction = str2double(inputMatrix{1}(impsat + 1)); %If 1 it uses the 2nd derivative of dfdS
            Sequential.ImplicitSolver.tol = str2double(inputMatrix{1}(impsat + 2));
            Sequential.ImplicitSolver.maxNewton = str2double(inputMatrix{1}(impsat + 3));
            %%%Initialise objects for outputting statistics
            Sequential.ImplicitSolver.timestep = 0;
            Sequential.ImplicitSolver.Newtons = 0;
            Sequential.ImplicitSolver.Chops = 0;
        end
        FIM = 0;
    case ('FIM')
        %%%%FIM settings
        FIM.MaxIter = str2double(inputMatrix{1}(settings + 1));
        FIM.Tol = str2double(inputMatrix{1}(settings + 2));
        FIM.CFL = str2double(inputMatrix{1}(settings + 3));
        %%%Initialise objects for outputting statistics
        FIM.timestep = zeros(TimeSteps,1);
        FIM.Iter = zeros(TimeSteps,1);
        FIM.Chops = zeros(TimeSteps,1);
        Sequential = 0;
end
if (str2double(inputMatrix{1}(adm + 1))~=0)
    ADMSettings.active = 1;
    temp = strfind(inputMatrix{1}, 'PRESSURE_INTERPOLATOR');
    x = find(~cellfun('isempty', temp));
    ADMSettings.Pressure_Interpolator =  char(inputMatrix{1}(x+1));
    temp = strfind(inputMatrix{1}, 'LEVELS');
    x = find(~cellfun('isempty', temp));
    ADMSettings.maxLevel = str2double(inputMatrix{1}(x+1));
    temp = strfind(inputMatrix{1}, 'TOLERANCE');
    x = find(~cellfun('isempty', temp));
    ADMSettings.tol = str2double(inputMatrix{1}(x+1));
    temp = strfind(inputMatrix{1}, 'COARSENING_RATIOS');
    x = find(~cellfun('isempty', temp));
    cx = str2double(inputMatrix{1}(x+1));
    cy = str2double(inputMatrix{1}(x+2));
    ADMSettings.Coarsening = [cx, cy; cx^2, cy^2; cx^3, cy^3]; %Coarsening Factors: Cx1, Cy1; Cx2, Cy2; ...; Cxn, Cyn;
else
    ADMSettings.active = 0;
end
end