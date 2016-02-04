%Define Wells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Inj, Prod] = WellsProperties(inj, prod, inputMatrix, Grid, K)
%Create objects 
Inj = struct('cells', {}, 'p', {}, 'r', {}, 'PI', {});
Prod = struct('cells', {}, 'p', {}, 'r', {}, 'PI', {});
%Injection wells
for i=1:length(inj)
    Inj(i).cells = Addwellcells(Grid.Nx, char(inputMatrix(inj(i) + 1)), str2double(inputMatrix(inj(i) + 2)), str2double(inputMatrix(inj(i) + 3)), str2double(inputMatrix(inj(i) + 4)));
    Inj(i).p = str2double(inputMatrix(inj(i) + 5)); %[Pa]
    Inj(i).r = str2double(inputMatrix(inj(i) + 6)); %Well radius in m
    Inj(i).PI = 2000;
    %Inj(i).PI = ComputeProductivityIndex(Inj.r, K(1, Inj(i).x, Inj(i).y), K(2, Inj(i).x, Inj(i).y), Grid.dx, Grid.dy, 1);
end
%Production Wells
for i=1:length(prod)
    Prod(i).cells = Addwellcells(Grid.Nx, char(inputMatrix(prod(i) + 1)), str2double(inputMatrix(prod(i) + 2)), str2double(inputMatrix(prod(i) + 3)), str2double(inputMatrix(prod(i) + 4)));
    Prod(i).p = str2double(inputMatrix(prod(i) + 5)); %[Pa]
    Prod(i).r = str2double(inputMatrix(prod(i) + 6)); %Well radius in m
    Prod(i).PI = 2000;
    %Prod(i).PI = ComputeProductivityIndex(Prod.r,K(1, Prod.x, Prod.y), K(2, Prod.x, Prod.y), Grid.dx, Grid.dy, 1);
end
end