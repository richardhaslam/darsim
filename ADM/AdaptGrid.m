%ADM - Grid construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DLGRGrid, ActiveFine, CoarseGrid, FineGrid] = AdaptGrid(FineGrid, CoarseGrid, S, Ro, Rw, maxLevel, tol)
% Coarsen the grid where resolution is not necessary
global NoWellsCoarseCells

%1. Choose Active Coarse cells and Flag fine ones
RefinementCriterion = 'SaturationDelta';

CoarseGrid(1).Active = NoWellsCoarseCells;
[FineGrid, CoarseGrid(1)] = ...
    SelectCoarseFine(FineGrid, CoarseGrid(1), 1, RefinementCriterion, S, Rw, Ro, tol);
for x = 2:maxLevel
    CoarseGrid(x).Active = DefinePossibleActive(CoarseGrid(x), CoarseGrid(x-1), x);
    [CoarseGrid(x-1), CoarseGrid(x)] = ...
        SelectCoarseFine(CoarseGrid(x-1), CoarseGrid(x), x, RefinementCriterion, S, Rw, Ro, tol);
end

%3. Count total number of active nodes
TotalActive = sum(FineGrid.Active);
NumberOfActive = zeros(maxLevel +1, 1);
NumberOfActive(1) = TotalActive;
for x = 1:maxLevel
    NumberOfActive(x+1) = sum(CoarseGrid(x).Active);
    TotalActive = TotalActive + sum(CoarseGrid(x).Active);
end
ActiveFine = reshape(FineGrid.Active, FineGrid.Nx, FineGrid.Ny);

%4. Add cells to DLGR grid
Field1 = 'N';  Value1 = NumberOfActive;
Field2 = 'I';  Value2 = zeros(TotalActive, 1);
Field3 = 'J';  Value3 = zeros(TotalActive, 1);
Field4 = 'CoarseFactor'; Value4 = zeros(TotalActive, 1); 
Field5 = 'CellIndex'; Value5 = zeros(TotalActive, 1);
Field6 = 'Father'; Value6 = zeros(TotalActive, maxLevel);
Field7 = 'Centers'; Value7 = zeros(TotalActive, maxLevel);
DLGRGrid = struct(Field1,Value1,Field2,Value2,Field3,Value3,Field4,Value4,Field5,Value5, Field6, Value6, Field7, Value7);

%Add fine Grid cells
[DLGRGrid, Tot] = AddActiveCells(DLGRGrid, FineGrid, 0, 0);

%Add Coarse Grids cells
for x = 1:maxLevel
    [DLGRGrid, Tot] = AddActiveCells(DLGRGrid, CoarseGrid(x), Tot, x);
end
DLGRGrid.maxLevel = max(DLGRGrid.level);
end