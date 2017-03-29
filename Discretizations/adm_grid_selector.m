%  ADM Grid Selector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector < handle
    properties
        tol
        NoWellsCoarseCells
    end
    methods
        function obj = adm_grid_selector(tol)
            obj.tol = tol;
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, maxLevel)
            % Select the ADM grid
            FineGrid.Active = ones(FineGrid.N, 1);
            
            % Coarsen the grid where resolution is not necessary
            S = reshape(ProductionSystem.Reservoir.State.z(:,1), FineGrid.Nx, FineGrid.Ny, FineGrid.Nz); % it's useful in this form for selecting the grid.
            
            %1. Choose Active Coarse cells and Flag fine ones  
            CoarseGrid(1).Active = obj.NoWellsCoarseCells;
            obj.SelectCoarseFine(FineGrid, CoarseGrid(1), 1, S);
            for x = 2:maxLevel
                obj.DefinePossibleActive(CoarseGrid(x), CoarseGrid(x-1), x);
                obj.SelectCoarseFine(CoarseGrid(x-1), CoarseGrid(x), x, S);
            end
            
            %3. Count total number of active nodes
            TotalActive = sum(FineGrid.Active);
            NumberOfActive = zeros(maxLevel +1, 1);
            NumberOfActive(1) = TotalActive;
            for x = 1:maxLevel
                NumberOfActive(x+1) = sum(CoarseGrid(x).Active);
                TotalActive = TotalActive + sum(CoarseGrid(x).Active);
            end
            
            ADMGrid.Initialize(TotalActive, NumberOfActive, maxLevel);
            
            %Add fine Grid cells
            obj.AddActiveCells(ADMGrid, FineGrid, 0);
            
            %Add Coarse Grids cells
            for x = 1:maxLevel
                obj.AddActiveCells(ADMGrid, CoarseGrid(x), x);
            end
            % I need to know the maximum coarse level that was used.
            ADMGrid.MaxLevel = max(ADMGrid.level);
        end
        function DefinePossibleActive(obj, CoarseGrid, FineGrid, level)
            % For a given level defines possible active cells
            
            CoarseGrid.Active = ones(CoarseGrid.N, 1);
            % If a cell inside the block is refined the whole block cannot be coarsened
            Nf = FineGrid.N;
            for i=1:Nf
                if FineGrid.Active(i) == 0
                    CoarseGrid.Active(FineGrid.Fathers(i, level)) = 0;
                    %Active(FineGrid.Father(FineGrid.Neighbours(i).indexes, level)) = 0;
                end
            end
            
            % Force the jump between two neighbouring cells to be max 1 level!
            Nc = CoarseGrid.N;
            temp = 1 - CoarseGrid.Active;
            for j=1:Nc
                if CoarseGrid.Active(j) == 1
                    vecNei = CoarseGrid.Neighbours(j).indexes;
                    check = sum(temp(vecNei));
                    if check > 0
                        CoarseGrid.Active(j) = 0;
                    end
                end
            end
        end
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, level, S)
            %Given a Fine and a Coarse Grids chooses the cells that have to be active
            %1. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                I = CoarseGrid.I(c,2);
                J = CoarseGrid.J(c,2);
                K = CoarseGrid.K(c,2);
                Imin = I - floor((CoarseGrid.CoarseFactor(1) - 1)/2);
                Imax = I + ceil((CoarseGrid.CoarseFactor(1) - 1)/2);
                Jmin = J - floor((CoarseGrid.CoarseFactor(2) - 1)/2);
                Jmax = J + ceil((CoarseGrid.CoarseFactor(2) - 1)/2);
                Kmin = K - floor((CoarseGrid.CoarseFactor(3) - 1)/2);
                Kmax = K + ceil((CoarseGrid.CoarseFactor(3) - 1)/2);
                
                % Max e Min saturation
                Smax = max(max(max(S(Imin:Imax, Jmin:Jmax, Kmin:Kmax))));
                Smin = min(min(min(S(Imin:Imax, Jmin:Jmax, Kmin:Kmax))));
                if CoarseGrid.Active(c) == 1
                    n = CoarseGrid.Neighbours(c).indexes;
                    Nn = length(n);
                    i = 1;
                    while i <= Nn
                        if (abs(Smax-S(CoarseGrid.I(n(i), 2), CoarseGrid.J(n(i), 2), CoarseGrid.K(n(i), 2)))...
                                > obj.tol || abs(Smin-S(CoarseGrid.I(n(i),2),CoarseGrid.J(n(i),2), CoarseGrid.K(n(i), 2))) > obj.tol)
                            CoarseGrid.Active(c) = 0;
                            %CoarseGrid.Active(i) = 0;
                            i = Nn + 1;
                        else
                            i = i+1;
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
            %Cindeces = find();
            FineGrid.Active(CoarseGrid.Children(CoarseGrid.Active == 1,:)) = 0;
        end
        function AddActiveCells(obj, ADMGrid, Grid, level)
            count = 0;
            for i=1:Grid.N
                if(Grid.Active(i) == 1)
                    h = ADMGrid.Ntot + count + 1;
                    ADMGrid.I(h) = Grid.I(i);
                    ADMGrid.J(h) = Grid.J(i);
                    ADMGrid.K(h) = Grid.K(i);
                    ADMGrid.CoarseFactor(h, 1) = Grid.CoarseFactor(1);
                    ADMGrid.CoarseFactor(h, 2) = Grid.CoarseFactor(2);
                    ADMGrid.CoarseFactor(h, 3) = Grid.CoarseFactor(3);
                    ADMGrid.CellIndex(h) = i;
                    ADMGrid.level(h) = level;
                    ADMGrid.Fathers(h, :) = Grid.Fathers(i, :);
                    ADMGrid.Children{h} = Grid.Children(i,:);
                    ADMGrid.GrandChildren{h} = Grid.GrandChildren(i,:);
                    ADMGrid.Verteces(h,:) = Grid.Verteces(i,:);
                    count = count + 1;
                end
            end
            ADMGrid.Ntot = ADMGrid.Ntot + count;
        end
    end
end