%  ADM Grid Selector: delta prop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 30 June 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector_temperature < adm_grid_selector
    properties
    end
    methods
        function obj = adm_grid_selector_temperature(tol, key)
            obj@adm_grid_selector(tol);
            obj.key = key;
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, Residual, maxLevel)
            % SELECT the ADM GRID for next time-step based on dela x
            %% 1. Reset all cells to be active and stor x{m} 
            n_media = length(FineGrid);
            Tf = cell(n_media, 1);
            for m=1:n_media
                FineGrid(m).Active = ones(FineGrid(m).N, 1);
                if m==1
                    % for now wells are only in the reservoir
                    CoarseGrid(m,1).Active = obj.NoWellsCoarseCells;
                    Tf{m} = ProductionSystem.Reservoir.State.Properties(obj.key).Value;
                else
                    CoarseGrid(m,1).Active = ones(CoarseGrid(m).N, 1);
                    Tf{m} = ProductionSystem.FracturesNetwork.Fractures(m-1).State.Properties(obj.key).Value;
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
                        obj.SelectCoarseFine(FineGrid(m), CoarseGrid(m, 1), Tf{m});
                    elseif l <= maxLevel(m)
                        obj.SelectCoarseFine(CoarseGrid(m, l-1), CoarseGrid(m, l), Tf{m});
                    else
                        CoarseGrid(m, l).Active = zeros(CoarseGrid(m, l).N, 1);
                    end
                end
            end
            obj.CreateADMGrid(ADMGrid, FineGrid, CoarseGrid, maxLevel);
        end
        function SelectCoarseFine(obj, FineGrid, CoarseGrid, T)
            %Given a Fine and a Coarse Grids chooses the cells that have to be active
            %1. Select Active Coarse Blocks
            Nc = CoarseGrid.N;
            for c = 1:Nc
                % temperature of fine_cells belonging to coarse block c
                Tf_children = T(CoarseGrid.GrandChildren{c, :});
                % Max e Min temperature inside c
                Tfmax = max(Tf_children);
                Tfmin = min(Tf_children);
                if CoarseGrid.Active(c) == 1
                    n = CoarseGrid.Neighbours(c).indexes;
                    Nn = length(n);
                    i = 1;
                    while i <= Nn
                        % Max min temperature of neighbour n(i)
                        Tf_children = T(CoarseGrid.GrandChildren{n(i), :});
                        Tfn_max = max(Tf_children);
                        Tfn_min = min(Tf_children);
                        if (abs(Tfmax-Tfn_min) > obj.tol || abs(Tfmin-Tfn_max) > obj.tol)
                            CoarseGrid.Active(c) = 0;
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
            FineGrid.Active([CoarseGrid.Children{CoarseGrid.Active == 1,:}]) = 0;
        end
    end
end