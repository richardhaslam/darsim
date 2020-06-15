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
        isCoupled
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
                %CoarseGrid(m).Active = ones(CoarseGrid(m).N, 1);
                CoarseGrid(m).Active(:) = 1;
                
                % If a cell inside a CoarseGrid block is refined, the whole block cannot be coarsened
                CoarseGrid(m).Active( unique( FineGrid(m).Fathers(FineGrid(m).Active == 0,1) ) ) = 0;
                
                % Force the jump between two neighbouring cells to be max 1 level!
                Nc = CoarseGrid(m).N;
                temp = 1 - CoarseGrid(m).Active;
                for j=1:Nc
                    if CoarseGrid(m).Active(j) == 1
                        vecNei = CoarseGrid(m).Neighbours{j};
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
                NumberOfActive(m, 1) = sum(FineGrid(m).Active);
                TotalActive = TotalActive + sum(FineGrid(m).Active);
                Nc_global(m, 1) = sum([FineGrid(1:m-1).N]);
                for x = 1:max(maxLevel)
                    Nc_global(m, x+1) = sum([CoarseGrid(1:m-1, x).N].*[CoarseGrid(1:m-1, x).hasCoarseNodes]);
                    if CoarseGrid(m,x).hasCoarseNodes
                        NumberOfActive(m, x+1) = sum(CoarseGrid(m, x).Active);
                        TotalActive = TotalActive + sum(CoarseGrid(m, x).Active);
                    end
                end
            end
            
            % 2. Initialise Grid
            ADMGrid.Initialize(TotalActive, NumberOfActive, max(maxLevel));
            
            % 3. Add fine Grid cells
            for m=1:n_media
                obj.AddActiveCells(ADMGrid, FineGrid(m), 0, Nc_global(m, :));
            end
            
            % 4. Add Coarse Grids cells
            for L = 1:max(maxLevel)
                for m=1:n_media
                    if CoarseGrid(m,L).hasCoarseNodes
                        obj.AddActiveCells(ADMGrid, CoarseGrid(m, L), L, Nc_global(m, :));
                    end
                end
            end
            % I need to know the maximum coarse level that was used.
            ADMGrid.MaxLevel = max(ADMGrid.level);
        end
        function AddActiveCells(obj, ADMGrid, Grid, level, N_global)
            count = 0;
            MaxStaticLevel = length(N_global)-1;
            for i=1:Grid.N
                if(Grid.Active(i) == 1)
                    h = ADMGrid.Ntot + count + 1;
                    ADMGrid.CoarseFactor(h, :) = Grid.CoarseFactor;
                    ADMGrid.CellIndex(h) = i + N_global(level+1); % The variable "N_global" is used to convert the CoarseGrid numbering into ADM Global numbering.
                    ADMGrid.level(h) = level;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if level == 0
                        % 1. Adding Fathers
                        if ~isempty(Grid.Fathers)
                            ADMGrid.Fathers(h, :) = Grid.Fathers(i, :) + N_global(2:end);
                        else
                            ADMGrid.Fathers(h, :) = -1; % Assigning an invalid dummy value
                        end
                        
                        % 2. Adding Children
                        % There is no Children in the finescale level as there is no level below it.
                        
                        % 3. Adding Verteces
                        if ~isempty(Grid.Verteces)
                            ADMGrid.Verteces(h,:) = Grid.Verteces(i,:);
                        else
                            ADMGrid.Verteces(h, :) = -1; % Assigning an invalid dummy value
                        end
                        
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    elseif level > 0  &&  level < MaxStaticLevel
                        % 1. Adding Fathers
                        if ~isempty(Grid.Fathers)
                            ADMGrid.Fathers(h, level+1:end) = Grid.Fathers(i, :) + N_global(level+2:end);
                            ADMGrid.Fathers(h, 1:level) = -1; % Assigning an invalid dummy value
                        else
                            ADMGrid.Fathers(h, :) = -1; % Assigning an invalid dummy value
                        end
                        
                        % 2. Adding Children
                        for L = 1 : level
                            ADMGrid.Children{h, L} = Grid.Children{i, L} + N_global(level-L+1);
                        end
                        
                        % 3. Adding Verteces
                        if ~isempty(Grid.Verteces)
                            ADMGrid.Verteces(h, level+1:end) = Grid.Verteces(i, :);
                            ADMGrid.Verteces(h, 1:level) = -1; % Assigning an invalid dummy value
                        else
                            ADMGrid.Verteces(h, :) = -1; % Assigning an invalid dummy value
                        end
                        
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    elseif level == MaxStaticLevel
                        % 1. Adding Fathers
                        % There is no Father in the coarsest level as there is no level above it.
                        ADMGrid.Fathers(h, :) = -1; % Assigning an invalid dummy value
                        
                        % 2. Adding Children
                        for L = 1 : level
                            ADMGrid.Children{h, L} = Grid.Children{i, L} + N_global(level-L+1);
                        end
                        
                        % 3. Adding Verteces
                        % There is no vertex in the coarsest level as there is no level above it.
                        ADMGrid.Verteces(h, :) = -1; % Assigning an invalid dummy value
                    end
                    
                    count = count + 1;
                end
            end
            ADMGrid.Ntot = ADMGrid.Ntot + count;
        end
        function AddActiveFineCells(obj, ADMGrid, Grid, level, N_global)
            ADMGrid.Fathers(h, :) = Grid.Fathers(i, :) + N_global(2:end);
        end
        function AddActiveCoaeseCells(obj, ADMGrid, Grid, level, N_global)
        end
    end
end