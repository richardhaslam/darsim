%ADM - Define possible active cells for level l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 1 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Active = DefinePossibleActive(CoarseGrid, FineGrid, level)
% For a given level defines possible active cells

Active = ones(CoarseGrid.Nx*CoarseGrid.Ny,1);
%If a cell inside the block is refined the whole block cannot be coarsened
Nf = FineGrid.Nx*FineGrid.Ny;
for i=1:Nf
    if FineGrid.Active(i) == 0
        Active(FineGrid.Father(i, level)) = 0;
        %Active(FineGrid.Father(FineGrid.Neighbours(i).indexes, level)) = 0;
    end
end

%Force the jump between two neighbouring cells to be max 1 level!
Nc = CoarseGrid.Nx*CoarseGrid.Ny;
temp = 1 - Active;
for j=1:Nc
    if Active(j) == 1
        vecNei = CoarseGrid.Neighbours(j).indexes;
        check = sum(temp(vecNei));
        if check > 0
            Active(j) = 0;
        end
    end
end
end