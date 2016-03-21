%ADM SET-UP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 21 March 2016
%Last modified: 21 March 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%Set-up all objects needed by ADM%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% Construct Coarse Grids
CoarseGridding;

%% Pressure interpolators
disp('Pressure interpolator - start computation');
CoarseGrid = PressureInterpolator(Grid, K, CoarseGrid, maxLevel, ADMSettings.Pressure_Interpolator);
disp('Pressure interpolator - end')
disp(char(2));

