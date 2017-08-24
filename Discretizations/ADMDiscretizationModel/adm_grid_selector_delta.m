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
            obj@adm_grid_selector(tol, key);
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, maxLevel)
            % SELECT the ADM GRID for next time-step
            FineGrid.Active = ones(FineGrid.N, 1);
            
            % Coarsen the grid where resolution is not necessary
            S = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
            
            % Choose Active Coarse cells and Flag fine ones  
            CoarseGrid(1).Active = obj.NoWellsCoarseCells;
            obj.SelectCoarseFine(FineGrid, CoarseGrid(1), S);
            for x = 2:maxLevel
                obj.DefinePossibleActive(CoarseGrid(x), CoarseGrid(x-1), x);
                obj.SelectCoarseFine(CoarseGrid(x-1), CoarseGrid(x), S);
            end
            
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, S)
            %Given a Fine and a Coarse Grids chooses the cells that have to be active
            %1. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                % Saturation of fine_cells belonging to coarse block c
                S_children = S(CoarseGrid.GrandChildren(c, :));
                % Max e Min saturation inside c
                Smax = max(S_children);
                Smin = min(S_children);
                if CoarseGrid.Active(c) == 1
                    n = CoarseGrid.Neighbours(c).indexes;
                    Nn = length(n);
                    i = 1;
                    while i <= Nn
                        % Man mix saturation of neighbour n(i)
                        S_children = S(CoarseGrid.GrandChildren(n(i), :));
                        Sn_max = max(S_children);
                        Sn_min = min(S_children);
                        if (abs(Smax-Sn_min) > obj.tol || abs(Smin-Sn_max) > obj.tol)
                            CoarseGrid.Active(c) = 0;
                            %CoarseGrid.Active(i) = 0;
                            i = Nn + 1;
                        else
                            i = i+1;
                        end
                    end
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