%Define Wells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Inj, Prod] = WellsProperties(inj, prod, inputMatrix, Grid, K)
%Wells

%Injection wells
Inj.r = str2double(inputMatrix(inj + 4)); %Well radius in m
Inj.p = str2double(inputMatrix(inj + 3)); %[Pa]
Inj.x = str2double(inputMatrix(inj + 1));
Inj.y = str2double(inputMatrix(inj + 2));
Inj.PI = ComputeProductivityIndex(Inj.r, K(1, Inj.x, Inj.y), K(2, Inj.x, Inj.y), Grid.dx, Grid.dy, 1);
%Production Wells
Prod.r = str2double(inputMatrix(prod + 4)); %Well radius in m
Prod.p = str2double(inputMatrix(prod + 3)); %[Pa]
Prod.x = str2double(inputMatrix(prod + 1));
Prod.y = str2double(inputMatrix(prod + 2));
Prod.PI = ComputeProductivityIndex(Prod.r, K(1, Prod.x, Prod.y), K(2, Prod.x, Prod.y), Grid.dx, Grid.dy, 1);
end