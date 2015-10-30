%Construct a coarse grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoarseGrid] = BuildCoarseGrid(Grid, CoarseGrid)
%Construct a coarse Grid given coarsening ration and fine Grid
CoarseGrid.Nx = Grid.Nx/CoarseGrid.CoarseFactor(1);
CoarseGrid.Ny = Grid.Ny/CoarseGrid.CoarseFactor(2);
Nc=CoarseGrid.Nx*CoarseGrid.Ny;
%Coordinates of the centres
CoarseGrid.I = ones(Nc, 1);
CoarseGrid.J = ones(Nc, 1);
Jindexes = ceil(CoarseGrid.CoarseFactor(2)/2):CoarseGrid.CoarseFactor(2):Grid.Ny;
for i=1:CoarseGrid.Ny
    a = CoarseGrid.Nx*(i-1)+1;
    CoarseGrid.I(a:a+CoarseGrid.Nx-1) = ceil(CoarseGrid.CoarseFactor(1)/2):CoarseGrid.CoarseFactor(1):Grid.Nx;
    CoarseGrid.J(a:a+CoarseGrid.Nx-1) = Jindexes(i)*ones(CoarseGrid.Nx,1);
end
CoarseGrid.Active = zeros(CoarseGrid.Nx*CoarseGrid.Ny,1);
CoarseGrid.Wells = zeros(CoarseGrid.Nx*CoarseGrid.Ny,1);
CoarseGrid = AssignNeighbours(CoarseGrid);
end