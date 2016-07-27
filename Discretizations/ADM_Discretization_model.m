%  ADM discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef ADM_Discretization_model < Discretization_model
    properties
        ADMSettings
        CoarseGrid
        ADMGrid
        NoWellsCoarseCells
    end
    methods
        function obj = ADM_Discretization_model(nx, ny, nz, settings)
            obj@Discretization_model(nx, ny, nz);
            obj.ADMSettings = settings;
            obj.CoarseGrid = coarse_grid();
        end
        function Initialize(obj, ProductionSystem)
            disp('Algebraic Dynamic Multilevel (ADM) method run')
            obj.ReservoirGrid.Initialize(ProductionSystem.Reservoir);
            % Perforated cells
            obj.DefinePerforatedCells(ProductionSystem.Wells);
            
            % Construct Coarse Grids
            obj.ConstructCoarseGrids(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            %% Pressure interpolators
            disp('Pressure interpolator - start computation');
            obj.StaticPressureInterpolators(ProductionSystem.Reservoir.K);
            disp('Pressure interpolator - end')
            disp(char(2));
        end
        function ConstructCoarseGrids(obj, Inj, Prod)
            % Add coordinates to fine-scale grid (useful for ADM)
            obj.ReservoirGrid.AddCoordinates();
            
            % Construct all coarse grids
            for i=1:obj.ADMSettings.maxLevel
                obj.CoarseGrid(i).CoarseFactor = obj.ADMSettings.Coarsening(i,:);
                obj.CoarseGrid(i).BuildCoarseGrid(obj.ReservoirGrid);
            end
            
            % Assign fathers and verteces 
            obj.AssignFathers();
            
            % Flag coarse blocks with wells
            obj.CoarseWells(Inj, Prod);
            obj.NoWellsCoarseCells = ones(obj.CoarseGrid(1).N, 1);
            Nc1 = obj.CoarseGrid(1).N;
            if obj.ADMSettings.maxLevel > 1
               for i = 1:Nc1
                   if obj.CoarseGrid(2).Wells(obj.CoarseGrid(1).Father(i,2)) == 1
                       obj.NoWellsCoarseCells(i) = 0;
                   end
               end
            else
               for i =1:Nc1
                   if obj.CoarseGrid(1).Wells(i) == 1
                       obj.NoWellsCoarseCells(i) = 0;
                   end
               end
            end
        end
        function AssignFathers(obj)
            % Fine Grid
            obj.ReservoirGrid.Father = zeros(obj.ReservoirGrid.N, obj.ADMSettings.maxLevel);
            obj.ReservoirGrid.Centers = zeros(obj.ReservoirGrid.N, obj.ADMSettings.maxLevel);
            for i=1:obj.ADMSettings.maxLevel
                Nc = obj.CoarseGrid(i).Nx*obj.CoarseGrid(i).Ny;
                for c = 1:Nc
                    Imin = obj.CoarseGrid(i).I(c) - floor((obj.CoarseGrid(i).CoarseFactor(1) - 1)/2);
                    Imax = obj.CoarseGrid(i).I(c) + ceil((obj.CoarseGrid(i).CoarseFactor(1) - 1)/2);
                    Jmin = obj.CoarseGrid(i).J(c) - floor((obj.CoarseGrid(i).CoarseFactor(2) - 1)/2);
                    Jmax = obj.CoarseGrid(i).J(c) + ceil((obj.CoarseGrid(i).CoarseFactor(2) - 1)/2);
                    fc = obj.CoarseGrid(i).I(c) + (obj.CoarseGrid(i).J(c) - 1)*obj.ReservoirGrid.Nx;
                    obj.ReservoirGrid.Centers(fc,i) = 1;
                    %Scan fine cells inside c
                    for I = Imin:Imax
                        for J = Jmin:Jmax
                            f = I + (J-1)*obj.ReservoirGrid.Nx;
                            obj.ReservoirGrid.Father(f,i) = c;
                        end
                    end
                end
                obj.CoarseGrid(i).Centers = ones(obj.CoarseGrid(i).Nx*obj.CoarseGrid(i).Ny, 1);
            end
            % Coarse Grids
            for x = 1:obj.ADMSettings.maxLevel-1
                Nf = obj.CoarseGrid(x).Nx*obj.CoarseGrid(x).Ny;
                obj.CoarseGrid(x).Father = zeros(Nf, obj.ADMSettings.maxLevel);
                obj.CoarseGrid(x).Centers = zeros(Nf, obj.ADMSettings.maxLevel);
                %obj.CoarseGrid(x).Centers(:, x) = ones(Nf, 1);
                for y = x+1:obj.ADMSettings.maxLevel
                    Nc = obj.CoarseGrid(y).Nx*obj.CoarseGrid(y).Ny;
                    for c = 1:Nc
                        Imin = obj.CoarseGrid(y).I(c) - floor((obj.CoarseGrid(y).CoarseFactor(1) - 1)/2);
                        Imax = obj.CoarseGrid(y).I(c) + ceil((obj.CoarseGrid(y).CoarseFactor(1) - 1)/2);
                        Jmin = obj.CoarseGrid(y).J(c) - floor((obj.CoarseGrid(y).CoarseFactor(2) - 1)/2);
                        Jmax = obj.CoarseGrid(y).J(c) + ceil((obj.CoarseGrid(y).CoarseFactor(2) - 1)/2);
                        for l = 1:Nf
                            if ((obj.CoarseGrid(x).I(l) <= Imax && obj.CoarseGrid(x).I(l) >= Imin) && ...
                                    (obj.CoarseGrid(x).J(l) <= Jmax && obj.CoarseGrid(x).J(l) >= Jmin))
                                obj.CoarseGrid(x).Father(l,y) = c;
                            end
                            if (obj.CoarseGrid(x).I(l) == obj.CoarseGrid(y).I(c) && obj.CoarseGrid(x).J(l) == obj.CoarseGrid(y).J(c))
                                obj.CoarseGrid(x).Centers(l, y) = 1;
                            end
                        end
                    end
                end
            end
            obj.CoarseGrid(obj.ADMSettings.maxLevel).Father = zeros(obj.CoarseGrid(obj.ADMSettings.maxLevel).N, obj.ADMSettings.maxLevel);
            obj.CoarseGrid(obj.ADMSettings.maxLevel).Centers = zeros(obj.CoarseGrid(obj.ADMSettings.maxLevel).N, obj.ADMSettings.maxLevel);
        end
        function CoarseWells(obj, Inj, Prod)
            for i=1:length(Inj)
                % Flag coarse Nodes with wells
                I = Inj(i).Cells;
                for x = 1:obj.ADMSettings.maxLevel
                    obj.CoarseGrid(x).Wells(obj.ReservoirGrid.Father(I,x)) = 1;
                end
            end
            for i =1:length(Prod)
                P = Prod(i).Cells;
                for x = 1:obj.ADMSettings.maxLevel
                    obj.CoarseGrid(x).Wells(obj.ReservoirGrid.Father(P,x)) = 1;
                end
            end
        end
        function StaticPressureInterpolators(obj, K)
            % Restriction and Prolongation for Pressure blocks
            switch (obj.ADMSettings.Pressure_Interpolator)
                case ('Constant')
                    %If it's constant interpolation I only need to build restriction
                    obj.CoarseGrid(1).MsR = MsRestriction(FineGrid, obj.CoarseGrid(1), FineGrid.Nx*FineGrid.Ny, obj.CoarseGrid(1).N, 1);
                    obj.CoarseGrid(1).MsP = obj.CoarseGrid(1).MsR';
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
                    [Grid] = ComputeTransmissibility(FineGrid, K);
                    Tx = Grid.Tx;
                    Ty = Grid.Ty;
                    %Assemble pressure matrix
                    Ap = AssemblePressureMatrix(Grid);
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
        
        function MsR = MsRestriction(FineGrid, CoarseGrid, Nf, Nc, level)
            %% MSFV Restriction Operator
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
    end
end