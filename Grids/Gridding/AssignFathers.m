% Assign fathers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 13 June 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Grid, CoarseGrid] = AssignFathers(Grid, CoarseGrid, maxLevel)
% Fine Grid
Grid.Father = zeros(Grid.N, maxLevel);
Grid.Centers = zeros(Grid.N, maxLevel);
for i=1:maxLevel
    Nc = CoarseGrid(i).Nx*CoarseGrid(i).Ny;
    for c = 1:Nc
        Imin = CoarseGrid(i).I(c) - floor((CoarseGrid(i).CoarseFactor(1) - 1)/2);
        Imax = CoarseGrid(i).I(c) + ceil((CoarseGrid(i).CoarseFactor(1) - 1)/2);
        Jmin = CoarseGrid(i).J(c) - floor((CoarseGrid(i).CoarseFactor(2) - 1)/2);
        Jmax = CoarseGrid(i).J(c) + ceil((CoarseGrid(i).CoarseFactor(2) - 1)/2);
        fc = CoarseGrid(i).I(c) + (CoarseGrid(i).J(c) - 1)*Grid.Nx;
        Grid.Centers(fc,i) = 1;
        %Scan fine cells inside c
        for I = Imin:Imax
            for J = Jmin:Jmax
                f = I + (J-1)*Grid.Nx;
                Grid.Father(f,i) = c;
            end
        end
    end
    CoarseGrid(i).Centers = ones(CoarseGrid(i).Nx*CoarseGrid(i).Ny, 1);
end
% Coarse Grids
for x = 1:maxLevel-1
    Nf = CoarseGrid(x).Nx*CoarseGrid(x).Ny;
    CoarseGrid(x).Father = zeros(Nf, maxLevel);
    CoarseGrid(x).Centers = zeros(Nf, maxLevel);
    %CoarseGrid(x).Centers(:, x) = ones(Nf, 1);
    for y = x+1:maxLevel
        Nc = CoarseGrid(y).Nx*CoarseGrid(y).Ny;
        for c = 1:Nc
            Imin = CoarseGrid(y).I(c) - floor((CoarseGrid(y).CoarseFactor(1) - 1)/2);
            Imax = CoarseGrid(y).I(c) + ceil((CoarseGrid(y).CoarseFactor(1) - 1)/2);
            Jmin = CoarseGrid(y).J(c) - floor((CoarseGrid(y).CoarseFactor(2) - 1)/2);
            Jmax = CoarseGrid(y).J(c) + ceil((CoarseGrid(y).CoarseFactor(2) - 1)/2);
            for l = 1:Nf
                if ((CoarseGrid(x).I(l) <= Imax && CoarseGrid(x).I(l) >= Imin) && ...
                        (CoarseGrid(x).J(l) <= Jmax && CoarseGrid(x).J(l) >= Jmin))
                    CoarseGrid(x).Father(l,y) = c;
                end
                if (CoarseGrid(x).I(l) == CoarseGrid(y).I(c) && CoarseGrid(x).J(l) == CoarseGrid(y).J(c))
                    CoarseGrid(x).Centers(l, y) = 1;
                end
            end
        end
    end
end
CoarseGrid(maxLevel).Father = zeros(CoarseGrid(maxLevel).Nx*CoarseGrid(maxLevel).Ny, maxLevel);
CoarseGrid(maxLevel).Centers = zeros(CoarseGrid(maxLevel).Nx*CoarseGrid(maxLevel).Ny, maxLevel);
end
%CoarseGrid(maxLevel).Centers(:, maxLevel) = ones(CoarseGrid(maxLevel).Nx*CoarseGrid(maxLevel).Ny, 1);