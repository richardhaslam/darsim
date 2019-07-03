%  ADM Grid Selector: delta prop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_delta2 < adm_grid_selector
    properties
        tol_time
    end
    methods
        function obj = adm_grid_selector_delta2(tol_space, tol_time, key)
            obj@adm_grid_selector(tol_space);
            obj.key = key;
            obj.tol_time = tol_time;
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, Residual, maxLevel)
            % SELECT the ADM GRID for next time-step based on delta x
            %% 1. Reset all cells to be active and stor x{m}
            n_media = length(FineGrid);
            delta = cell(n_media, 1);
            S = cell(n_media, 1);
            for m=1:n_media
                FineGrid(m).Active = ones(FineGrid(m).N, 1);
                if m==1
                    % for now wells are only in the reservoir
                    CoarseGrid(m,1).Active = obj.NoWellsCoarseCells;
                    num = ProductionSystem.Reservoir.State.Properties(obj.key).Value - ...
                        ProductionSystem.Reservoir.State_old.Properties(obj.key).Value;
                    delta{m} = num;
                    S{m} = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
                    
                else
                    CoarseGrid(m,1).Active = ones(CoarseGrid(m).N, 1);
                    num = ProductionSystem.FracturesNetwork.Fractures(m-1).State.Properties(obj.key).Value - ...
                        ProductionSystem.FracturesNetwork.Fractures(m-1).State_old.Properties(obj.key).Value;
                    delta{m} = num;
                    S{m} = ProductionSystem.FracturesNetwork.Fractures(m-1).State.Properties(obj.key).Value;
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
                        obj.SelectCoarseFine(FineGrid(m), CoarseGrid(m, 1), delta{m}, S{m},  false);
                    elseif l <= maxLevel(m)
                        obj.SelectCoarseFine(CoarseGrid(m, l-1), CoarseGrid(m, l), delta{m}, S{m}, true);
                    else
                        CoarseGrid(m, l).Active = zeros(CoarseGrid(m, l).N, 1);
                    end
                end
            end
             %% 3. Create ADM Grid
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        
        
        function SelectGridCoarse(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, Residual, maxLevel)
            % SELECT the ADM GRID for next time-step based on delta x
            %% 1. Reset all cells to be active and stor x{m}
            n_media = length(FineGrid);
            for m=1:n_media
                FineGrid(m).Active = ones(FineGrid(m).N, 1);
                if m==1
                    % for now wells are only in the reservoir
                    CoarseGrid(m,1).Active = obj.NoWellsCoarseCells;
                else
                    CoarseGrid(m,1).Active = ones(CoarseGrid(m).N, 1);
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
                        obj.SelectCoarseFineC(FineGrid(m), CoarseGrid(m, 1));
                    elseif l <= maxLevel(m)
                        obj.SelectCoarseFineC(CoarseGrid(m, l-1), CoarseGrid(m, l));
                    else
                        CoarseGrid(m, l).Active = zeros(CoarseGrid(m, l).N, 1);
                    end
                end
            end
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, delta, S, flag)
            %Given a Fine and a Coarse Grid chooses the cells that have to be active
            %1. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                % Saturation of fine_cells belonging to coarse block c
                S_children = S([CoarseGrid.GrandChildren{c, :}]);
                % Max e Min saturation inside c
                Smax = max(S_children);
                Smin = min(S_children);
                if CoarseGrid.Active(c) == 1
                    n = CoarseGrid.Neighbours(c).indexes;
                    Nn = length(n);
                    i = 1;
                    while i <= Nn
                        % Man mix saturation of neighbour n(i)
                        S_children = S([CoarseGrid.GrandChildren{n(i), :}]);
                        Sn_max = max(S_children);
                        Sn_min = min(S_children);
                        if flag == false
                            if (abs(Smax-Sn_min) > obj.tol || abs(Smin-Sn_max) > obj.tol)
                                CoarseGrid.Active(c) = 0;
                                %CoarseGrid.Active(i) = 0 
                                if max(delta([CoarseGrid.GrandChildren{c, :}]))< obj.tol_time
                                    CoarseGrid.Active(c) = 1;
                                end
                                i = Nn + 1;
                            else
                                i = i+1;
                            end
                        else
                            if (abs(Smax-Sn_min) > 2*obj.tol || abs(Smin-Sn_max) > 2*obj.tol)
                                CoarseGrid.Active(c) = 0;
                                %CoarseGrid.Active(i) = 0;
                                i = Nn + 1;
                            else
                                i = i+1;
                            end
                        end
                    end
                end
            end
            
%             %2. Do not coarsen neighbours of cells that are fine
%                         DummyActive = CoarseGrid.Active;
%                         for i = 1:Nc
%                             if (CoarseGrid.Active(i) == 0)
%                                 DummyActive(CoarseGrid.Neighbours(i).indexes) = 0;
%                             end
%                         end
%                         CoarseGrid.Active = DummyActive.*CoarseGrid.Active;
            FineGrid.Active([CoarseGrid.Children{CoarseGrid.Active == 1,:}]) = 0;
        end
        function SelectCoarseFineC(obj, FineGrid, CoarseGrid)
            FineGrid.Active([CoarseGrid.Children{CoarseGrid.Active == 1,:}]) = 0;
        end
    end
end