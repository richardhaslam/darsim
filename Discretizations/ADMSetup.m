%ADM SET-UP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 21 March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%Set-up all objects needed by ADM%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Grid, CoarseGrid] = ADMSetup(Grid, K, ADMSettings, Inj, Prod)
%% Construct Coarse Grids
Grid.CoarseFactor = [1, 1];
[Grid, CoarseGrid] = CoarseGridding(Grid, ADMSettings, Inj, Prod);

%% Pressure interpolators
disp('Pressure interpolator - start computation');
CoarseGrid = PressureInterpolator(Grid, K, CoarseGrid, ADMSettings.maxLevel, ADMSettings.Pressure_Interpolator);
disp('Pressure interpolator - end')
disp(char(2));
end
