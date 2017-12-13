%  ADM Grid Selector: time criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 June 2017
%Last modified: 27 November 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_time < adm_grid_selector
    properties
        tol2 = .2;
        Epsilon_old;
    end
    methods
        function obj = adm_grid_selector_time(tol, key, N, maxLevels)
            obj@adm_grid_selector(tol);
            obj.key = key;
            obj.Epsilon_old = zeros(N, maxLevels);
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, Residual, maxLevel)
            % SELECT the ADM GRID for next time-step
            % Grid is chosen based on (deltaX)^{n} = X^{n} - X^{n-1}
            
            %% 1. Reset all cells to be active and compute (deltaX)^{n} 
            n_media = length(FineGrid);
            delta = cell(n_media, 1);
            S = cell(n_media, 2);
            for m=1:n_media
                FineGrid(m).Active = ones(FineGrid(m).N, 1);
                if m==1
                    % for now wells are only in the reservoir
                    CoarseGrid(m,1).Active = obj.NoWellsCoarseCells;
                    num = ProductionSystem.Reservoir.State.Properties(obj.key).Value - ...
                        ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
                    delta{m} = num;
                    S{m, 1} = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
                    S{m, 2} = ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
                else
                    CoarseGrid(m,1).Active = ones(CoarseGrid(m).N, 1);
                    num = ProductionSystem.FracturesNetwork.Fractures(m-1).State.Properties(obj.key).Value - ...
                        ProductionSystem.FracturesNetwork.Fractures(m-1).State_old.Properties(obj.key).Value;
                    delta{m} = num;
                    S{m, 1} = ProductionSystem.FracturesNetwork.Fractures(m-1).State.Properties(obj.key).Value;
                    S{m, 2} = ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
                end
            end
                            
            %% 2. Select active cells
            Tol = obj.tol;
            Tol2 = obj.tol2;
            for l=1:max(maxLevel)
                % 2.a choose possible active grids for level l
                if l>1 
                    obj.DefinePossibleActive(CoarseGrid(:, l), CoarseGrid(:, l-1), l);
                end
                for m=1:n_media
                    % 2.b choose active cells of level l 
                    if l==1
                        % coarse grid 1 to 0 (fine-scale)
                        obj.SelectCoarseFine(FineGrid(m), CoarseGrid(m, 1), delta{m}, S{m,:}, l);
                    elseif l <= maxLevel(m)
                        obj.SelectCoarseFine(CoarseGrid(m, l-1), CoarseGrid(m, l), delta{m}, S{m,:}, l);
                    else
                        CoarseGrid(m, l).Active = zeros(CoarseGrid(m, l).N, 1);
                    end
                end
                obj.tol = obj.tol / 2;
                %obj.tol2 = obj.tol2 / 2;
            end
            obj.tol = Tol;
            obj.tol2 = Tol2;
            
            %% 3. Create ADM Grid
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, delta, S, S_old, l)
            %Given a Fine (level l-1) and a Coarse (level l) Grids chooses the cells that have to be active
            
            %% 2. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                % fine-scale cells inside coarse block c
                indexes_fs = CoarseGrid.GrandChildren(c,:);        
 
                Max = max(delta(indexes_fs)); Max(abs(Max)<1e-3) = 0;
                Min = min(delta(indexes_fs)); Min(abs(Min)<1e-3) = 1;
                DeltaS = norm(delta(indexes_fs), inf);
                Epsilon = (S(indexes_fs) - S_old(indexes_fs)) ./ (mean(S(indexes_fs)) - mean(S_old(indexes_fs)));
                Gradient = norm(Epsilon - obj.Epsilon_old(indexes_fs, l), inf);
                if CoarseGrid.Active(c) == 1 && Gradient > 1e-3 && sum(~isnan(Epsilon)) && norm(S_old(indexes_fs), inf) > 0.101 
                    CoarseGrid.Active(c) = 0;
                elseif CoarseGrid.Active(c) == 1 && ( (CoarseGrid.DeltaS(c) < 1e-3 && DeltaS > 1e-3) || Max*Min<0)
                    CoarseGrid.Active(c) = 0;
                end
                CoarseGrid.DeltaS(c) = DeltaS;
                obj.Epsilon_old(indexes_fs, l) = (S(indexes_fs) - S_old(indexes_fs)) ./ (mean(S(indexes_fs)) - mean(S_old(indexes_fs)));
            end
            
            %% 3. Set to inactive fine blocks (level l-1) belonging to active Coarse Blocks (level l)
            FineGrid.Active(CoarseGrid.Children(CoarseGrid.Active == 1,:)) = 0;
        end
    end
end