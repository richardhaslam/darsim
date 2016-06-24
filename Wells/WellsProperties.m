%Define Wells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Inj, Prod] = WellsProperties(inj, prod, inputMatrix, Grid, K)
%Injection wells
for i=1:length(inj)
    Well = Addwellcells(Grid.Nx, char(inputMatrix(inj(i) + 1)), str2double(inputMatrix(inj(i) + 2)), str2double(inputMatrix(inj(i) + 3)), str2double(inputMatrix(inj(i) + 4)));
    if strcmp(char(inputMatrix(inj(i) + 5)), 'rate')
        Well.type = 'RateConstrained';
        Well.q = str2double(inputMatrix(inj(i) + 6))*Grid.Lx*Grid.Ly*Grid.h*Grid.por/(24*3600);
    else
        Well.type = 'PressureConstrained';
        Well.p = str2double(inputMatrix(inj(i) + 6)); %[Pa]
    end
    Well.r = str2double(inputMatrix(inj(i) + 6)); %Well radius in m
    Well.PI = 2000;
    Inj(i) = Well;
    Inj(i).z = [0, 0];
    Inj(i).s = 0;
    Inj(i).x1 = [0, 0];
    %Inj(i).PI = ComputeProductivityIndex(Inj.r, K(1, 1, 1), K(2, 1, 1), Grid.dx, Grid.dy, 1);
end
%Production Wells
for i=1:length(prod)
    Well = Addwellcells(Grid.Nx, char(inputMatrix(prod(i) + 1)), str2double(inputMatrix(prod(i) + 2)), str2double(inputMatrix(prod(i) + 3)), str2double(inputMatrix(prod(i) + 4)));
    if strcmp(char(inputMatrix(prod(i) + 5)), 'rate')
        Well.type = 'RateConstrained';
        Well.q = - str2double(inputMatrix(prod(i) + 6))*Grid.Lx*Grid.Ly*Grid.h*Grid.por/(24*3600);
    else
        Well.type = 'PressureConstrained';
        Well.p = str2double(inputMatrix(prod(i) + 6)); %[Pa]
    end
    Well.r = str2double(inputMatrix(prod(i) + 6)); %Well radius in m
    Well.PI = 2000;
    Prod(i) = Well;
    %Prod(i).PI = ComputeProductivityIndex(Prod.r,K(1, 216, 54), K(2, 216, 54), Grid.dx, Grid.dy, 1);
end
end
%% Perforated cells
function [Well] = Addwellcells(Nx, Direction, I, Y, final)
global Errors
switch Direction
    case ('vertical')
        Well.cells = I + (Y-1)*Nx;
    case ('horizontal')
        Yindexes = Y:1:final;
        Well.cells = ones(1, length(Yindexes))*I + (Yindexes - 1)*Nx;
    otherwise
        disp('ERROR: unknown well direction!! Please select either vertical or horizontal');
        Errors = 1;
        Well = 'error';
end
end