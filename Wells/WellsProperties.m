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
Inj = addwellcoordinates(str2double(inputMatrix(inj + 1)), str2double(inputMatrix(inj + 2)), str2double(inputMatrix(inj + 3)), str2double(inputMatrix(inj + 4)));
Inj.p = str2double(inputMatrix(inj + 5)); %[Pa]
Inj.r = str2double(inputMatrix(inj + 6)); %Well radius in m
%Inj.PI = ComputeProductivityIndex(Inj.r, K(1, Inj.x, Inj.y), K(2, Inj.x, Inj.y), Grid.dx, Grid.dy, 1);
Inj.PI = 2;
%Production Wells
Prod = addwellcoordinates(str2double(inputMatrix(prod + 1)), str2double(inputMatrix(prod + 2)), str2double(inputMatrix(prod + 3)), str2double(inputMatrix(prod + 4)));
Prod.p = str2double(inputMatrix(prod + 5)); %[Pa]
Prod.r = str2double(inputMatrix(prod + 6)); %Well radius in m
%Prod.PI = ComputeProductivityIndex(Prod.r,K(1, Prod.x, Prod.y), K(2, Prod.x, Prod.y), Grid.dx, Grid.dy, 1);
Prod.PI = 2;
end