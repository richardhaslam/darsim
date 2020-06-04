%  ADM Grid Selector: delta prop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 30 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_delta < adm_grid_selector
    properties
    end
    methods
        function obj = adm_grid_selector_delta(tol, key)
            obj@adm_grid_selector(tol);
            obj.key = key;
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, Residual, maxLevel)
            % SELECT the ADM GRID for next time-step based on delta x
            %% 1. Reset all cells to be active and store Var{m}
            n_media = length(FineGrid);
            Var = cell(n_media, 1);
            for m=1:n_media
                FineGrid(m).Active(:) = 1;
                if m==1
                    % for now wells are only in the reservoir
                    CoarseGrid(m,1).Active = obj.NoWellsCoarseCells;
                    Var{m} = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
                else
                    CoarseGrid(m,1).Active(:) = 1;
                    Var{m} = ProductionSystem.FracturesNetwork.Fractures(m-1).State.Properties(obj.key).Value;
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
                        obj.SelectCoarseFine(FineGrid(m), CoarseGrid(m, 1), Var{m}, false);
                    elseif l <= maxLevel(m)
                        obj.SelectCoarseFine(CoarseGrid(m, l-1), CoarseGrid(m, l), Var{m}, true);
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
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, Var, flag)
            %Given a Fine and a Coarse Grid chooses the cells that have to be active
            %1. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                % Values of the variable "Var" for FineGrid cells belonging to CoarseGrid block c
                Var_children = Var(CoarseGrid.Children{c, end});
                % Finding max and min values inside CoarseGrid block c
                VarMax = max(Var_children);
                VarMin = min(Var_children);
                if CoarseGrid.Active(c) == 1
                    Neighbor = CoarseGrid.Neighbours(c).indexes;
                    NofNeighbor = length(Neighbor);
                    i = 1;
                    while i <= NofNeighbor
                        % Finding max and min values of neighbour i
                        Var_children_Neighbor = Var(CoarseGrid.Children{Neighbor(i), end});
                        VarNeighborMax = max(Var_children_Neighbor);
                        VarNeighborMin = min(Var_children_Neighbor);
                        if flag == false
                            if (abs(VarMax-VarNeighborMin) > obj.tol || abs(VarMin-VarNeighborMax) > obj.tol)
                                CoarseGrid.Active(c) = 0;
                                i = NofNeighbor + 1;
                            else
                                i = i+1;
                            end
                        else
                            if (abs(VarMax-VarNeighborMin) > 2*obj.tol || abs(VarMin-VarNeighborMax) > 2*obj.tol)
                                CoarseGrid.Active(c) = 0;
                                i = NofNeighbor + 1;
                            else
                                i = i+1;
                            end
                        end
                    end
                end
            end
            %2. Do not coarsen neighbours of cells that are at finescale resolution
            DummyActive = CoarseGrid.Active;
            for i = 1:Nc
                if (CoarseGrid.Active(i) == 0)
                    DummyActive(CoarseGrid.Neighbours(i).indexes) = 0;
                end
            end
            CoarseGrid.Active = DummyActive.*CoarseGrid.Active;
            
            %3. Set to inactive the fine block belonging to Active Coarse Blocks
            %if CoarseGrid.N ~= 0
                FineGrid.Active([CoarseGrid.Children{CoarseGrid.Active == 1,1}]) = 0;
            %end
        end
        function SelectCoarseFineC(obj, FineGrid, CoarseGrid)
            FineGrid.Active([CoarseGrid.Children{CoarseGrid.Active == 1,:}]) = 0;
        end
    end
end