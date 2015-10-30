%ADM - Define possible active cells for level l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Active = DefinePossibleActive(CoarseGrid, FineGrid, level)
% For a given level defines possible active cells
Active = ones(CoarseGrid.Nx*CoarseGrid.Ny,1);
Nf = FineGrid.Nx*FineGrid.Ny;
for i=1:Nf
    if FineGrid.Active(i) == 0
        Active(FineGrid.Father(i, level)) = 0;
    end
end
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