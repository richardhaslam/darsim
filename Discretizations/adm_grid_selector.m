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
        function CreateADMGrid(obj, ADMGrid, FineGrid, CoarseGrid, maxLevel)
            % 1. Count total number of active nodes
            TotalActive = sum(FineGrid.Active);
            NumberOfActive = zeros(maxLevel +1, 1);
            NumberOfActive(1) = TotalActive;
            for x = 1:maxLevel
                NumberOfActive(x+1) = sum(CoarseGrid(x).Active);
                TotalActive = TotalActive + sum(CoarseGrid(x).Active);
            end
            
            % 2. Initialise Grid
            ADMGrid.Initialize(TotalActive, NumberOfActive, maxLevel);
            
            % 3. Add fine Grid cells
            obj.AddActiveCells(ADMGrid, FineGrid, 0);
            
            % 4. Add Coarse Grids cells
            for x = 1:maxLevel
                obj.AddActiveCells(ADMGrid, CoarseGrid(x), x);
            end
            % I need to know the maximum coarse level that was used.
            ADMGrid.MaxLevel = max(ADMGrid.level);
        end
        function AddActiveCells(obj, ADMGrid, Grid, level)
            count = 0;
            for i=1:Grid.N
                if(Grid.Active(i) == 1)
                    h = ADMGrid.Ntot + count + 1;
                    ADMGrid.I(h) = Grid.I(i);
                    ADMGrid.J(h) = Grid.J(i);
                    ADMGrid.K(h) = Grid.K(i);
                    ADMGrid.CoarseFactor(h, 1) = Grid.CoarseFactor(1);
                    ADMGrid.CoarseFactor(h, 2) = Grid.CoarseFactor(2);
                    ADMGrid.CoarseFactor(h, 3) = Grid.CoarseFactor(3);
                    ADMGrid.CellIndex(h) = i;
                    ADMGrid.level(h) = level;
                    ADMGrid.Fathers(h, :) = Grid.Fathers(i, :);
                    ADMGrid.Children{h} = Grid.Children(i,:);
                    ADMGrid.GrandChildren{h} = Grid.GrandChildren(i,:);
                    ADMGrid.Verteces(h,:) = Grid.Verteces(i,:);
                    count = count + 1;
                end
            end
            ADMGrid.Ntot = ADMGrid.Ntot + count;
        end
    end
end