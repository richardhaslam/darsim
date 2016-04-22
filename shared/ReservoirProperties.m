%Define Reservoir Properties, Grid and wells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Grid, K] = ReservoirProperties(size, grid, perm, pert, por, inputMatrix)
%Dimensions
Grid.Lx = str2double(inputMatrix(size +1));  %Dimension in x−direction [m] 
Grid.Ly = str2double(inputMatrix(size +2));  %Dimension in y−direction [m]
Grid.h = str2double(inputMatrix(size + 3));  %Reservoir thickness [m]

%Gridding
Grid.Nx = str2double(inputMatrix(grid +1));  %Nx
Grid.dx = Grid.Lx/Grid.Nx;                   %dx
Grid.Ny = str2double(inputMatrix(grid +2));  %Ny
Grid.dy = Grid.Ly/Grid.Ny;                   %dy
Grid.Ax = Grid.dy*Grid.h;                    %Cross section in x direction
Grid.Ay = Grid.dx*Grid.h;                    %Cross section in y direction
Grid.Volume = Grid.dx.*Grid.dy*Grid.h;       %Cell volume [m^3]
Grid.por = str2double(inputMatrix(por + 1)); %Porosity
Grid.N = Grid.Nx*Grid.Ny;                    %Total number of cells  

%Rock permeability in [m^2].
if strcmp(inputMatrix(perm - 1), 'INCLUDE')
    file  = strcat('../Permeability/', char(inputMatrix(perm +1))); %File name
    field = load(file);     %load the file in a vector
    field = reshape(field(3:end),[field(1) field(2)]);     % reshape it to specified size
    Kx = reshape(field(1:Grid.Nx,1:Grid.Ny)*10^(-12), Grid.N, 1);  % make it the size of the grid
    Ky = Kx;
    K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
elseif (~isempty(pert))
    %Perturbation parameters
    Grid.Vdp = str2double(inputMatrix(pert +2));
    Grid.Lcx = str2double(inputMatrix(pert +4));
    Grid.Lcy = str2double(inputMatrix(pert +6));
    Grid.Kav = str2double(inputMatrix(pert +8));
    
    [K] = Randomfield(Grid);
    Kx = K;
    Ky = K; 
    Kx = reshape(Kx,Grid.N,1);
    Ky = reshape(Ky,Grid.N,1);
    K  = reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);    
else
    Kx = ones(Grid.N,1)*10^(-12);
    Ky = ones(Grid.N,1)*10^(-12);
    K=reshape([Kx, Ky]', 2, Grid.Nx, Grid.Ny);
end
end