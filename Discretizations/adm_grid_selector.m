%  ADM Grid Selector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef adm_grid_selector < handle
    properties
        tol
        NoWellsCoarseCells
    end
    methods
        function obj = adm_grid_selector(tol)
            obj.tol = tol;
        end
        function SelectGrid(obj, FineGrid, CoarseGrid, ADMGrid, ProductionSystem, maxLevel)
            % Coarsen the grid where resolution is not necessary
            
            %1. Choose Active Coarse cells and Flag fine ones  
            CoarseGrid(1).Active = obj.NoWellsCoarseCells;
            obj.SelectCoarseFine(FineGrid, CoarseGrid(1), 1, RefinementCriterion, S, Rw, Ro, tol);
            for x = 2:maxLevel
                obj.DefinePossibleActive(CoarseGrid(x), CoarseGrid(x-1), x);
                obj.SelectCoarseFine(CoarseGrid(x-1), CoarseGrid(x), x, RefinementCriterion, S, Rw, Ro, tol);
            end
            
            %3. Count total number of active nodes
            TotalActive = sum(FineGrid.Active);
            NumberOfActive = zeros(maxLevel +1, 1);
            NumberOfActive(1) = TotalActive;
            for x = 1:maxLevel
                NumberOfActive(x+1) = sum(CoarseGrid(x).Active);
                TotalActive = TotalActive + sum(CoarseGrid(x).Active);
            end
            
            %4. Add cells to DLGR grid
            %Field1 = 'N';  Value1 = NumberOfActive;
            %Field2 = 'I';  Value2 = zeros(TotalActive, 1);
            %Field3 = 'J';  Value3 = zeros(TotalActive, 1);
            %Field4 = 'CoarseFactor'; Value4 = zeros(TotalActive, 1);
            %Field5 = 'CellIndex'; Value5 = zeros(TotalActive, 1);
            %Field6 = 'Father'; Value6 = zeros(TotalActive, maxLevel);
            %Field7 = 'Centers'; Value7 = zeros(TotalActive, maxLevel);
            %DLGRGrid = struct(Field1,Value1,Field2,Value2,Field3,Value3,Field4,Value4,Field5,Value5, Field6, Value6, Field7, Value7);
            ADMGrid.initialize(TotalAcive);
            
            %Add fine Grid cells
            obj.AddActiveCells(ADMGrid, FineGrid, 0, 0);
            
            %Add Coarse Grids cells
            for x = 1:maxLevel
                obj.AddActiveCells(ADMGrid, CoarseGrid(x), Tot, x);
            end
            % I need to know the maximum coarse level that was used.
        end
    end
end