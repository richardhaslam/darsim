%Assign Wells to coarse cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoarseGrid] = CoarseWells(FineGrid, CoarseGrid, maxLevel, Inj, Prod)
for i=1:length(Inj)
    % Flag coarse Nodes with wells
    I = Inj(i).cells;
    for x = 1:maxLevel
        CoarseGrid(x).Wells(FineGrid.Father(I,x)) = 1;
    end
end
for i =1:length(Prod)
    P = Prod(i).cells; 
    for x = 1:maxLevel
        CoarseGrid(x).Wells(FineGrid.Father(P,x)) = 1;
    end
end
end