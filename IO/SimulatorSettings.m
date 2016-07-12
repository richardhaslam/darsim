%SIMULATOR SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Coupling, CouplingStats, ADMSettings, FlashSettings] = ...
    SimulatorSettings(TimeSteps, Strategy, settings, impsat, adm, inputMatrix)

%SIMULATOR SETTINGS
switch (Strategy)
    case ('Sequential')
        %%%%Sequential strategy settings
        Sequential.MaxExtIter = str2double(inputMatrix{1}(settings + 1));
        Sequential.Tol = str2double(inputMatrix{1}(settings +2));
        Sequential.CFL = str2double(inputMatrix{1}(settings +3));  %CFL number
        Sequential.ImpSat = 0;
        if (~isempty(impsat))
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
    case ('FIM')
        %%%%FIM settings
        NLSolver = NL_Solver;
        Coupling = FIM_Strategy('FIM', NLSolver);
        Coupling.CFL = str2double(inputMatrix{1}(settings + 3));
        
        Coupling.NLSolver.MaxIter = str2double(inputMatrix{1}(settings + 1));
        Coupling.NLSolver.Tol = str2double(inputMatrix{1}(settings + 2));
        
        %%%Initialise objects for outputting statistics
        CouplingStats = FIM_Stats(TimeSteps);
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

temp = strfind(inputMatrix{1}, 'FLASH LOOPS');
FlashSet = find(~cellfun('isempty', temp));
FlashSettings.TolInner=str2double(inputMatrix{1}(FlashSet + 2));
FlashSettings.MaxIt=str2double(inputMatrix{1}(FlashSet + 3));
FlashSettings.TolFlash=str2double(inputMatrix{1}(FlashSet + 4));
end