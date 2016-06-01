%ADM - Select active blocks for a given level l
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FineGrid, CoarseGrid] = SelectCoarseFine(FineGrid, CoarseGrid, level, RefinementCriterion, S, Rw, Ro, tol)
%Given a Fine and a Coarse Grids chooses the cells that have to be active
%1. Select Active Coarse Blocks
switch (RefinementCriterion)
    case('SaturationDelta')
        Nc = CoarseGrid.Nx*CoarseGrid.Ny;
        for c = 1:Nc
            I = CoarseGrid.I(c);
            J = CoarseGrid.J(c);
            Imin = I - floor((CoarseGrid.CoarseFactor(1) - 1)/2);
            Imax = I + ceil((CoarseGrid.CoarseFactor(1) - 1)/2);
            Jmin = J - floor((CoarseGrid.CoarseFactor(2) - 1)/2);
            Jmax = J + ceil((CoarseGrid.CoarseFactor(2) - 1)/2);
            %Max e Min saturation
            Smax = max(max(S(Imin:Imax, Jmin:Jmax)));
            Smin = min(min(S(Imin:Imax, Jmin:Jmax)));
            if CoarseGrid.Active(c) == 1
                n = CoarseGrid.Neighbours(c).indexes;
                Nn = length(n);
                i = 1;
                while i <= Nn
                    if (abs(Smax-S(CoarseGrid.I(n(i)),CoarseGrid.J(n(i))))...
                            > tol || abs(Smin-S(CoarseGrid.I(n(i)),CoarseGrid.J(n(i)))) > tol)
                        CoarseGrid.Active(c) = 0;
                        %CoarseGrid.Active(i) = 0;
                        i = Nn + 1;
                    else
                        i = i+1;
                    end
                end
            end
        end
    case('Residual')
        Nc = CoarseGrid.Nx*CoarseGrid.Ny;
        Rw = abs(Rw)./max(abs(Rw));
        Ro = abs(Ro)./max(abs(Ro));
        tol = 0.1;
        for c = 1:Nc
            I = CoarseGrid.I(c);
            J = CoarseGrid.J(c);
            Imin = I - floor((CoarseGrid.CoarseFactor(1) - 1)/2);
            Imax = I + ceil((CoarseGrid.CoarseFactor(1) - 1)/2);
            Jmin = J - floor((CoarseGrid.CoarseFactor(2) - 1)/2);
            Jmax = J + ceil((CoarseGrid.CoarseFactor(2) - 1)/2);
            i=Imin:Imax;
            j=Jmin:Jmax;
            [p,q] = meshgrid(i, j);
            pairs = [p(:), q(:)];
            %indexes of the fine cells
            indexes = pairs(:,1) + (pairs(:,2)-1)*CoarseGrid.Nx*CoarseGrid.CoarseFactor(1);
            if (CoarseGrid.Wells(c) == 1)
                CoarseGrid.Active(c) = 0;
            elseif CoarseGrid.Active(c) == 1
                if sum(Ro(indexes)) > tol*10
                    CoarseGrid.Active(c) = 0;
                end
            end
        end         
end

%2. Do not coarsen neighbors of cells that are fine
DummyActive = CoarseGrid.Active;
for i = 1:Nc
    if (CoarseGrid.Active(i) == 0)
        DummyActive(CoarseGrid.Neighbours(i).indexes) = 0;
    end
end
CoarseGrid.Active = DummyActive.*CoarseGrid.Active;

%3. Set to inactive fine block belonging to Active Coarse Blocks
Nf = FineGrid.Nx*FineGrid.Ny;
for i = 1:Nf
    if (CoarseGrid.Active(FineGrid.Father(i, level)) == 1)
        FineGrid.Active(i) = 0;
    end
end