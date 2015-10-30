%Assign Wells to coarse cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoarseGrid] = CoarseWells(FineGrid, CoarseGrid, maxLevel, Inj, Prod)
% Flag coarse Nodes with wells
I = Inj.x + (Inj.y - 1)*FineGrid.Nx;
P = Prod.x + (Prod.y - 1)*FineGrid.Nx;
for x = 1:maxLevel
    CoarseGrid(x).Wells(FineGrid.Father(I,x)) = 1;
    CoarseGrid(x).Wells(FineGrid.Father(P,x)) = 1;
end