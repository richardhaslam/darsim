%ADM - Pressure Interpolator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matteo Cusini's Research Code
%Author: Matteo Cusini
%TU Delft
%Created: 2015
%Last modified: 27 May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CoarseGrid = PressureInterpolator(FineGrid, Kt, CoarseGrid, maxLevel, Pressure_Interpolator)
% Restriction and Prolongation for Pressure blocks
switch (Pressure_Interpolator)
    case ('Constant')
        %If it's constant interpolation I only need to build restriction
        CoarseGrid(1).MsR = MsRestriction(FineGrid, CoarseGrid(1), FineGrid.Nx*FineGrid.Ny, CoarseGrid(1).Nx*CoarseGrid(1).Ny, 1);
        CoarseGrid(1).MsP = CoarseGrid(1).MsR';
        for x = 2:maxLevel
            CoarseGrid(x).MsR = MsRestriction(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x-1).Nx*CoarseGrid(x-1).Ny, CoarseGrid(x).Nx*CoarseGrid(x).Ny, x);
            CoarseGrid(x).MsP = CoarseGrid(x).MsR'; 
            CoarseGrid(x).constant = 1;
        end
    case ('Homogeneous')
        CoarseGrid(1).constant = 0;
        %Average permeability
        K = zeros(2,FineGrid.Nx, FineGrid.Ny);
        K(1,:,:) = ones(FineGrid.Nx, FineGrid.Ny)*mean(mean(Kt(1,:,:)));
        K(2,:,:) = ones(FineGrid.Nx, FineGrid.Ny)*mean(mean(Kt(2,:,:)));
        [Tx, Ty] = ComputeTransmissibility(FineGrid, K);
        %Assemble pressure matrix
        Ap = AssemblePressureMatrix(Tx, Ty, FineGrid.Nx, FineGrid.Ny);
        %Build MS operators
        [CoarseGrid(1).MsR, CoarseGrid(1).MsP, CoarseGrid(1).C] = MSFVOperators(FineGrid, CoarseGrid(1), Ap, 1);
        %Build first coarse system (with MsFV)
        CoarseGrid(1).A_c = CoarseGrid(1).MsR*Ap*CoarseGrid(1).MsP;
        for x = 2:maxLevel
            %Build MS operators
            [CoarseGrid(x).MsR, CoarseGrid(x).MsP, CoarseGrid(x).C] = MSFVOperators(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x-1).A_c, x);
            %Build coarse system (with MsFV)
            CoarseGrid(x).A_c = CoarseGrid(x).MsR*CoarseGrid(x-1).A_c*CoarseGrid(x).MsP;
            CoarseGrid(x).constant = 0;
        end
    case ('MS')
        CoarseGrid(1).constant = 0;
        %Reduce contrast to avoid peaks
        Ktvector = reshape(Kt(1,:, :),FineGrid.N, 1);
        lambdaMax = max(Ktvector);
        PlotPermeability(Kt, FineGrid);
        Ktvector(Ktvector./lambdaMax < 10^-2) = 10^-2*lambdaMax;
        Kt(1, :, :) = reshape(Ktvector, FineGrid.Nx, FineGrid.Ny);
        Kt(2, :, :) = reshape(Ktvector, FineGrid.Nx, FineGrid.Ny);
        PlotPermeability(Kt, FineGrid);
        %Compute transmissibility
        [Tx, Ty] = ComputeTransmissibility(FineGrid, Kt);
        %Assemble pressure matrix
        Ap = AssemblePressureMatrix(Tx, Ty, FineGrid.Nx, FineGrid.Ny);
        %Build MS operators
        [CoarseGrid(1).MsR, CoarseGrid(1).MsP, CoarseGrid(1).C] = MSFVOperators(FineGrid, CoarseGrid(1), Ap, 1);
        %Build first coarse system (with MsFE)
        CoarseGrid(1).A_c = CoarseGrid(1).MsP'*Ap*CoarseGrid(1).MsP;
        for x = 2:maxLevel
            %Build MS opeartors
            [CoarseGrid(x).MsR, CoarseGrid(x).MsP, CoarseGrid(x).C] = MSFVOperators(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x-1).A_c, x);
            %Construct coarse system (with MsFE)
            CoarseGrid(x).A_c = CoarseGrid(x).MsP'*CoarseGrid(x-1).A_c*CoarseGrid(x).MsP;
            CoarseGrid(x).constant = 0;
        end
end
end

%% MSFV Restriction Operator
function MsR = MsRestriction(FineGrid, CoarseGrid, Nf, Nc, level)
MsR = zeros(Nc, Nf);
if (FineGrid.CoarseFactor(1) == 1)
    for r=1:Nc
        %coordinates of fine cells contained in the coarse block
        Imin = CoarseGrid.I(r) - floor((CoarseGrid.CoarseFactor(1) - 1)/2);
        Imax = CoarseGrid.I(r) + ceil((CoarseGrid.CoarseFactor(1) - 1)/2);
        Jmin = CoarseGrid.J(r) - floor((CoarseGrid.CoarseFactor(2) - 1)/2);
        Jmax = CoarseGrid.J(r) + ceil((CoarseGrid.CoarseFactor(2) - 1)/2);
        i=Imin:Imax;
        j=Jmin:Jmax;
        [p,q] = meshgrid(i, j);
        pairs = [p(:), q(:)];
        %indexes of the fine cells
        c = pairs(:,1) + (pairs(:,2)-1)*FineGrid.Nx;
        %I make 1 those columns
        MsR(r,c) = 1;
    end
else
    for r=1:Nc
        for c=1:Nf
            if (r == FineGrid.Father(c, level))
                %I make 1 those columns
                MsR(r,c) = 1;
            end
        end
    end
end
MsR = sparse(MsR);
end