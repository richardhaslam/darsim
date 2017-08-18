%  ADM Grid Selector base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 June 2017
%Last modified: 30 June 2017
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
        function obj = adm_grid_selector(tol, key)
            obj.tol = tol;
            obj.key = key;
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
            for m=1:n_media
                TotalActive = TotalActive + sum(FineGrid(m).Active);
                NumberOfActive(m, 1) = sum(FineGrid(m).Active);
                for x = 1:max(maxLevel)
                    NumberOfActive(m, x+1) = sum(CoarseGrid(m, x).Active);
                    TotalActive = TotalActive + sum(CoarseGrid(m, x).Active);
                end
            end
            
            % 2. Initialise Grid
            ADMGrid.Initialize(TotalActive, NumberOfActive, max(maxLevel));
            
            % 3. Add fine Grid cells
            for m=1:n_media
                Nf = sum([FineGrid(1:m-1).N]);
                Nc = Nf;
                obj.AddActiveCells(ADMGrid, FineGrid(m), 0, Nf, Nc);
            end
            % 4. Add Coarse Grids cells
            for l = 1:max(maxLevel)
                for m=1:n_media
                    Nf = sum([FineGrid(1:m-1).N]);
                    if l <= maxLevel(m)
                        Nc = sum([CoarseGrid(1:m-1, l).N]);
                        obj.AddActiveCells(ADMGrid, CoarseGrid(m, l), l, Nf, Nc);
                    end
                end
            end
            % I need to know the maximum coarse level that was used.
            ADMGrid.MaxLevel = max(ADMGrid.level);
        end
        function AddActiveCells(obj, ADMGrid, Grid, level, Nf, Nc)
            count = 0;
            for i=1:Grid.N
                if(Grid.Active(i) == 1)
                    h = ADMGrid.Ntot + count + 1;
                    ADMGrid.I(h) = Grid.I(i);
                    ADMGrid.J(h) = Grid.J(i);
                    ADMGrid.K(h) = Grid.K(i);
                    ADMGrid.CoarseFactor(h, :) = Grid.CoarseFactor;
                    ADMGrid.CellIndex(h) = i + Nc;
                    ADMGrid.level(h) = level;
                    ADMGrid.Fathers(h, :) = Grid.Fathers(i, :);
                    ADMGrid.Children{h} = Grid.Children(i,:);
                    ADMGrid.GrandChildren{h} = Grid.GrandChildren(i,:) + Nf;
                    ADMGrid.Verteces(h,:) = Grid.Verteces(i,:);
                    count = count + 1;
                end
            end
            ADMGrid.Ntot = ADMGrid.Ntot + count;
        end
    end
end