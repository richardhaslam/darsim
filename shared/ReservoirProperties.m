%Define Reservoir Properties, Grid and wells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Grid, K] = ReservoirProperties(size, grid, perm, inputMatrix)
%Dimensions
Grid.Lx = str2double(inputMatrix(size +1));                              %Dimension in x−direction [m] 
Grid.Ly = str2double(inputMatrix(size +2));                              %Dimension in y−direction [m]
h = str2double(inputMatrix(size + 3));                                   %Reservoir thickness [m]

%Gridding
Grid.Nx = str2double(inputMatrix(grid +1)); 
Grid.dx = Grid.Lx/Grid.Nx; 
Grid.Ny = str2double(inputMatrix(grid +2)); 
Grid.dy = Grid.Ly/Grid.Ny; 
Grid.Ax = Grid.dy*h;                    %Cross section in x direction
Grid.Ay = Grid.dx*h;                    %Cross section in y direction
Grid.Volume = Grid.dx.*Grid.dy*h;       %Cell volume [m^3]
Grid.por = 0.2;                         %Porosity
Grid.N = Grid.Nx*Grid.Ny;  

%Rock permeability in [m^2].
if strcmp(inputMatrix(perm - 1), 'INCLUDE')
    file  = strcat('../Permeability/', char(inputMatrix(perm +1)));      % HERE YOU MENTION THE NAME OF THE FILE
    field = load(file);     % HERE YOU LOAD IT
    fieldX = field(2:2:length(field));      % HEAR YOU JUST READ THE KX ONES
    fieldX = reshape(fieldX,[60 220])';     % HERE YOU MAKE IT 60x220
    Kx = reshape(fieldX(1:Grid.Nx,1:Grid.Ny)*10^(-12), Grid.N, 1);  % HERE YOU TAKE THE SIZE YOU LIKE
    Ky =  Kx;
    K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
else
    Kx = ones(Grid.N,1)*10^(-12);
    Ky = ones(Grid.N,1)*10^(-12);
    K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
end
end