%COARSE Grids construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global NoWellsCoarseCells
Grid.CoarseFactor = [1, 1];

% Add I, J coordinates to Grid
Grid.I = ones(Grid.N, 1);
Grid.J = ones(Grid.N, 1);
Jindexes = 1:1:Grid.Ny;
for i=1:Grid.Ny
    a = Grid.Nx*(i-1)+1;
    Grid.I(a:a+Grid.Nx-1) = 1:1:Grid.Nx;
    Grid.J(a:a+Grid.Nx-1) = Jindexes(i)*ones(Grid.Nx,1);
end

%Build Coarse Grids
maxLevel = ADMSettings.maxLevel;
CoarseGrid = struct('CoarseFactor', {}, 'Nx', {}, 'Ny', {}, ...
    'I', {}, 'J', {},'Father', {}, 'Active', {}, 'Wells', {}, 'Neighbours', {}, 'Centers', {});
for i=1:maxLevel
    CoarseGrid(i).CoarseFactor = ADMSettings.Coarsening(i,:);
    CoarseGrid(i) = BuildCoarseGrid(Grid, CoarseGrid(i));
end
AssignFathers;

%Flag coarse blocks with wells
[CoarseGrid] = CoarseWells(Grid, CoarseGrid, maxLevel, Inj, Prod);
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