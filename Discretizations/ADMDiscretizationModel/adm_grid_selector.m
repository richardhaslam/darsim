%  ADM Grid Selector base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 June 2017
%Last modified: 18 December 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector < handle
    properties
        tol
        NoWellsCoarseCells
        key
    end
    methods (Abstract)
        obj = SelectGrid(obj)
    end
    methods
        function obj = adm_grid_selector(tol)
            obj.tol = tol;
        end
        function Initialise(obj,ProductionSystem, FineGrid, n_phases)
            % virtual call
        end
        function DefinePossibleActive(obj, CoarseGrid, FineGrid, level)
            % For a given level defines possible active cells
            for m=1:length(FineGrid)
                CoarseGrid(m).Active = ones(CoarseGrid(m).N, 1);
                % If a cell inside the block is refined the whole block cannot be coarsened
                Nf = FineGrid(m).N;
                for i=1:Nf
                    if FineGrid(m).Active(i) == 0
                        CoarseGrid(m).Active(FineGrid(m).Fathers(i, level)) = 0;
                        %Active(FineGrid.Father(FineGrid.Neighbours(i).indexes, level)) = 0;
                    end
                end
                
                % Force the jump between two neighbouring cells to be max 1 level!
                Nc = CoarseGrid(m).N;
                temp = 1 - CoarseGrid(m).Active;
                for j=1:Nc
                    if CoarseGrid(m).Active(j) == 1
                        vecNei = CoarseGrid(m).Neighbours(j).indexes;
                        check = sum(temp(vecNei));
                        if check > 0
                            CoarseGrid(m).Active(j) = 0;
                        end
                    end
                end
            end
        end
        function CreateADMGrid(obj, ADMGrid, FineGrid, CoarseGrid, maxLevel)
            n_media = length(maxLevel);
            NumberOfActive = zeros(n_media, max(maxLevel)+1);
            TotalActive = 0;
            
            % 1. Count total number of active nodes
            Nc_global = zeros(n_media, max(maxLevel) + 1);
            for m=1:n_media
                TotalActive = TotalActive + sum(FineGrid(m).Active);
                NumberOfActive(m, 1) = sum(FineGrid(m).Active);
                Nc_global(m, 1) = sum([FineGrid(1:m-1).N]);
                for x = 1:max(maxLevel)
                    Nc_global(m, x+1) = sum([CoarseGrid(1:m-1, x).N]);
                    NumberOfActive(m, x+1) = sum(CoarseGrid(m, x).Active);
                    TotalActive = TotalActive + sum(CoarseGrid(m, x).Active);
                end
            end
            
            % 2. Initialise Grid
            ADMGrid.Initialize(TotalActive, NumberOfActive, max(maxLevel));
            
            % 3. Add fine Grid cells
            for m=1:n_media
                obj.AddActiveCells(ADMGrid, FineGrid(m), 0, Nc_global(m, :));
            end
            % 4. Add Coarse Grids cells
            for l = 1:max(maxLevel)
                for m=1:n_media
                    if l <= maxLevel(m)
                        obj.AddActiveCells(ADMGrid, CoarseGrid(m, l), l, Nc_global(m, :));
                    end
                end
            end
            % I need to know the maximum coarse level that was used.
            ADMGrid.MaxLevel = max(ADMGrid.level);
        end
        function AddActiveCells(obj, ADMGrid, Grid, level, N_global)
            count = 0;
            for i=1:Grid.N
                if(Grid.Active(i) == 1)
                    h = ADMGrid.Ntot + count + 1;
                    ADMGrid.CoarseFactor(h, :) = Grid.CoarseFactor;
                    ADMGrid.CellIndex(h) = i + N_global(level+1); % add current level global numbering
                    ADMGrid.level(h) = level;
                    ADMGrid.Fathers(h, :) = Grid.Fathers(i, :) + N_global(2:end); % add all coarse levels global numbering
                    ADMGrid.Children{h} = Grid.Children{i, :} + N_global(max(level, 1)); % add level l-1 global numbering
                    ADMGrid.GrandChildren{h} = Grid.GrandChildren{i, :} + N_global(1); % add fine-scale global numbering
                    ADMGrid.Verteces(h,:) = Grid.Verteces(i,:);
                    count = count + 1;
                end
            end
            ADMGrid.Ntot = ADMGrid.Ntot + count;
        end
    end
end