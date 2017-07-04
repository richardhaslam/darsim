%  ADM Grid Selector: time criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 June 2017
%Last modified: 4 July 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_time < adm_grid_selector
    properties
    end
    methods
        function obj = adm_grid_selector_time(tol, key)
            obj@adm_grid_selector(tol, key);
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, maxLevel)
            % SELECT the ADM GRID for next time-step
            % Grid is chosen based on (deltaX)^n
            
            %% 0. Reset all cells to be active
            FineGrid.Active = ones(FineGrid.N, 1);
            for i=1:maxLevel
                CoarseGrid(i).Active = ones(CoarseGrid(i).N, 1);
            end
            
            %% 1. Compute change of property X over previous time-step
            delta = abs(ProductionSystem.Reservoir.State.Properties(obj.key).Value - ...
                    ProductionSystem.Reservoir.State_old.Properties(obj.key).Value);
            delta = reshape(delta, FineGrid.Nx, FineGrid.Ny, FineGrid.Nz);    
            
            %% 2. Go from coarse to fine
            % 2.a coarse grids lmax to 2
            for i=maxLevel:-1:2
                obj.SelectCoarseFine(CoarseGrid(i-1), CoarseGrid(i), delta);
            end
            % 2.b coarse grid 1 to 0 (fine-scale)
            obj.SelectCoarseFine(FineGrid, CoarseGrid(1), delta);
            
            %% 3. Create ADM Grid
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid);
        end
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, delta)
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
                deltaMax = max(max(max(delta(Imin:Imax, Jmin:Jmax, Kmin:Kmax))));
                if CoarseGrid.Active(c) == 1 && deltaMax > obj.Tol
                   CoarseGrid.Active(c) = 0;
                end
            end
            
            %2. Do not coarsen neighbours of cells that are fine
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
    end
end