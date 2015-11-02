%COARSE Grids construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global NoWellsCoarseCells

% I copy the base grid to a new object
FineGrid = Grid;
%clear Grid;
FineGrid.CoarseFactor = [1, 1];

% Add I, J coordinates to FineGrid
FineGrid.I = ones(FineGrid.N, 1);
FineGrid.J = ones(FineGrid.N, 1);
Jindexes = 1:1:FineGrid.Ny;
for i=1:FineGrid.Ny
    a = FineGrid.Nx*(i-1)+1;
    FineGrid.I(a:a+FineGrid.Nx-1) = 1:1:FineGrid.Nx;
    FineGrid.J(a:a+FineGrid.Nx-1) = Jindexes(i)*ones(FineGrid.Nx,1);
end

%Build Coarse Grids
maxLevel = ADMSettings.maxLevel;
CoarseGrid = struct('CoarseFactor', {}, 'Nx', {}, 'Ny', {}, ...
    'I', {}, 'J', {},'Father', {}, 'Active', {}, 'Wells', {}, 'Neighbours', {}, 'Centers', {});
for i=1:maxLevel
    CoarseGrid(i).CoarseFactor = ADMSettings.Coarsening(i,:);
    CoarseGrid(i) = BuildCoarseGrid(FineGrid, CoarseGrid(i));
end
AssignFathers;

%Flag coarse blocks with wells
[CoarseGrid] = CoarseWells(FineGrid, CoarseGrid, maxLevel, Inj, Prod);
NoWellsCoarseCells = ones(CoarseGrid(1).Nx*CoarseGrid(1).Ny,1);
Nc1 = CoarseGrid(1).Nx * CoarseGrid(1).Ny;
if maxLevel > 1
for i = 1:Nc1
    if CoarseGrid(2).Wells(CoarseGrid(1).Father(i,2)) == 1
       NoWellsCoarseCells(i) = 0;
    end
end
else
    for i =1:Nc1
        if CoarseGrid(1).Wells(i) == 1
           NoWellsCoarseCells(i) = 0;
        end
    end
end