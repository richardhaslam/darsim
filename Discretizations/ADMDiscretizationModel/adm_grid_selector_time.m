%  ADM Grid Selector: time criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 June 2017
%Last modified: 21 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_time < adm_grid_selector
    properties
        tol2 = 5;
    end
    methods
        function obj = adm_grid_selector_time(tol, key)
            obj@adm_grid_selector(tol, key);
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, maxLevel)
            % SELECT the ADM GRID for next time-step
            % Grid is chosen based on (deltaX)^{n} = X^{n} - X^{n-1}
            
            %% 1. Reset all cells to be active and compute (deltaX)^{n} 
            n_media = length(FineGrid);
            delta = cell(n_media, 1);
            for m=1:n_media
                FineGrid(m).Active = ones(FineGrid(m).N, 1);
                if m==1
                    % for now wells are only in the reservoir
                    CoarseGrid(m,1).Active = obj.NoWellsCoarseCells;
                    num = ProductionSystem.Reservoir.State.Properties(obj.key).Value - ...
                        ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
                    delta{m} = num;
                else
                    CoarseGrid(m,1).Active = ones(CoarseGrid(m).N, 1);
                    num = ProductionSystem.FracturesNetwork.Fractures(m-1).State.Properties(obj.key).Value - ...
                        ProductionSystem.FracturesNetwork.Fractures(m-1).State_old.Properties(obj.key).Value;
                    delta{m} = num;
                end
            end
                            
            %% 2. Select active cells
            for l=1:max(maxLevel)
                % 2.a choose possible active grids for level l
                if l>1 
                    obj.DefinePossibleActive(CoarseGrid(:, l), CoarseGrid(:, l-1), l);
                end
                for m=1:n_media
                    % 2.b choose active cells of level l 
                    if l==1
                        % coarse grid 1 to 0 (fine-scale)
                        obj.SelectCoarseFine(FineGrid(m), CoarseGrid(m, 1), delta{m});
                    elseif l <= maxLevel(m)
                        obj.SelectCoarseFine(CoarseGrid(m, l-1), CoarseGrid(m, l), delta{m});
                    else
                        CoarseGrid(m, l).Active = zeros(CoarseGrid(m, l).N, 1);
                    end
                end
            end
            
            %% 3. Create ADM Grid
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, delta)
            %Given a Fine (level l-1) and a Coarse (level l) Grids chooses the cells that have to be active
            
            %% 2. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                % fine-scale cells inside coarse block c
                indexes_fs = CoarseGrid.GrandChildren(c,:);
                
                % Max delta inside block c
                deltaSum = sum(delta(indexes_fs));
                Max = max(delta(indexes_fs)); Max(abs(Max)<1e-4) = 0;
                Min = min(delta(indexes_fs)); Min(abs(Min)<1e-4) = 1; 
                deltaRatio = Max/Min;   
                CoarseGrid.DeltaS(c) = deltaSum;
                if CoarseGrid.Active(c) == 1 && (abs(deltaSum) > obj.tol || deltaRatio > obj.tol2 || deltaRatio<0)
                   CoarseGrid.Active(c) = 0;
                end
            end
            
            %% 3. Set to inactive fine blocks (level l-1) belonging to active Coarse Blocks (level l)
            FineGrid.Active(CoarseGrid.Children(CoarseGrid.Active == 1,:)) = 0;
        end
    end
end