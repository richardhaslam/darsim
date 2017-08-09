%  ADM Grid Selector: time criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 June 2017
%Last modified: 5 July 2017
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
            CoarseGrid(1).Active = obj.NoWellsCoarseCells;
            
            %% 1. Compute change of property X over previous time-step
            num = abs(ProductionSystem.Reservoir.State.Properties(obj.key).Value - ...
                  ProductionSystem.Reservoir.State_old.Properties(obj.key).Value);
            den = ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
            %den = 1;
            a = abs(den) < 1e-3;
            b = num > 1e-3;
            num(a+b == 2) = 1; 
            den(a == 1) = 1;
            delta = num ./ den;
                                       
            %% 2. Go from coarse to fine
            % 2.a coarse grid 1 to 0 (fine-scale)
            obj.SelectCoarseFine(FineGrid, CoarseGrid(1), delta);
            % 2.b coarse grids lmax to 2
            for i=2:maxLevel
                obj.DefinePossibleActive(CoarseGrid(i), CoarseGrid(i-1), i);
                obj.SelectCoarseFine(CoarseGrid(i-1), CoarseGrid(i), delta);
            end
            
            % 2.c keep perforated cells at the fine-scale
            %obj.RefineWells(ProductionSystem.Wells, FineGrid, CoarseGrid, maxLevel)
            
            %% 3. Create ADM Grid
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, delta)
            %Given a Fine (level l-1) and a Coarse (level l) Grids chooses the cells that have to be active
            
            %% 1. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                % fine-scale cells inside coarse block c
                indexes_fs = CoarseGrid.GrandChildren(c,:);
                
                % Max delta inside block c
                deltaSum = sum(delta(indexes_fs));
                deltaMax = max(delta(indexes_fs));
                CoarseGrid.DeltaS(c) = deltaSum;
                if CoarseGrid.Active(c) == 1 && deltaSum > obj.tol
                   CoarseGrid.Active(c) = 0;
                end
            end
%             dummy = CoarseGrid.Active;
%             for c=1:Nc
%                 if dummy(c) == 0
%                     CoarseGrid.Active(CoarseGrid.Neighbours(c).indexes) = 0;
%                 end
%             end
            %% 3. Set to inactive fine blocks (level l-1) belonging to active Coarse Blocks (level l)
            FineGrid.Active(CoarseGrid.Children(CoarseGrid.Active == 1,:)) = 0;
        end
        function RefineWells(obj, Wells, FineGrid, CoarseGrid, maxLevel)
%             for w = 1:length(Wells.Inj)       
%                     FineGrid.Active(Wells.Inj(w).Cells) = 1;
%             end
%             for w = 1:length(Wells.Prod)
%                 FineGrid.Active(Wells.Prod(w).Cells) = 1;
%             end
              indexes = CoarseGrid(1).Children(obj.NoWellsCoarseCells == 0, :);
              FineGrid.Active(indexes) = 1;
        end
    end
end