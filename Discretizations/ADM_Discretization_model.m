%  ADM discretization model arrayfun(@(x)coarse_grid(), 1:n, 'UniformOutput',false)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef ADM_Discretization_model < Discretization_model
    properties
        maxLevel
        Coarsening
        GridMapper
        CoarseGrid
        ADMGrid
        ADMGridSelector
        OperatorsHandler
        ADMStats
    end
    methods
        function obj = ADM_Discretization_model(maxlevel, coarsening)
            n = length(maxlevel);
            obj.CoarseGrid = arrayfun(@(x)coarse_grid(), 1:n, 'UniformOutput',false)';
            obj.ADMGrid    = arrayfun(@(x)adm_grid()   , 1:n, 'UniformOutput',false)';
            obj.GridMapper = grid_mapper();
            obj.Coarsening = coarsening;
            obj.maxLevel   = maxlevel;
        end
        function AddOperatorsHandler(obj, operatorshandler)
            obj.OperatorsHandler = operatorshandler;
        end
        function AddADMGridSelector(obj, gridselector)
            obj.ADMGridSelector = gridselector;
        end
        function InitializeMapping(obj, ProductionSystem, FluidModel)
            disp('Algebraic Dynamic Multilevel (ADM) method run')
            % Construct Coarse Grids
            disp(char(2));
            disp('Constructing coarse grids');
            obj.ConstructCoarseGrids(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            %% Pressure interpolators
            disp('Static operators - start computation');
            start = tic;
            obj.OperatorsHandler.ProlongationBuilders(1).BuildStaticOperators(ProductionSystem,...
                FluidModel, obj.ReservoirGrid, obj.FracturesGrid, obj.CrossConnections, obj.maxLevel, obj.CoarseGrid);
            disp('Static operators - end')
            timer = toc(start);
            disp(['Static operators construction took ', num2str(timer)])
            disp(char(2));
        end
        function ConstructCoarseGrids(obj, Inj, Prod)
            % Construct all coarse grids for reservoir
            obj.CoarseGrid{1}(1).CoarseFactor = obj.Coarsening{1}(1,:);
            obj.CoarseGrid{1}(1).BuildCoarseGrid(obj.ReservoirGrid);
            obj.GridMapper.BuildFamily(obj.CoarseGrid{1}(1), obj.ReservoirGrid, obj.Coarsening{1}(1,:), 1);
            for i=2:obj.maxLevel(1)
                obj.CoarseGrid{1}(i).CoarseFactor = obj.Coarsening{1}(i,:);
                obj.CoarseGrid{1}(i).BuildCoarseGrid(obj.ReservoirGrid);
                obj.GridMapper.BuildFamily(obj.CoarseGrid{1}(i), obj.CoarseGrid{1}(i-1), obj.Coarsening{1}(1,:), i);
            end
            
            % Fathers and Verteces for reservoir
            obj.GridMapper.AssignFathersandVerteces(obj.ReservoirGrid, obj.CoarseGrid{1}, obj.maxLevel(1))

            % Flag coarse blocks with wells for reservoir
            obj.CoarseWells(Inj, Prod);
            obj.ADMGridSelector.NoWellsCoarseCells = ones(obj.CoarseGrid{1}(1).N, 1);
            Nc1 = obj.CoarseGrid{1}(1).N;
            if obj.maxLevel(1) > 1
               for i = 1:Nc1
                   if obj.CoarseGrid{1}(2).Wells(obj.CoarseGrid{1}(1).Fathers(i, 2)) == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            else
               for i =1:Nc1
                   if obj.CoarseGrid{1}(1).Wells(i) == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            end
            
            % Construct all coarse grids for fractures
            for f = 1 : length(obj.maxLevel) - 1
                obj.CoarseGrid{1+f}(1).CoarseFactor = obj.Coarsening{1+f}(1,:);
                obj.CoarseGrid{1+f}(1).BuildCoarseGrid(obj.FracturesGrid.Grids(f));
                obj.GridMapper.BuildFamily(obj.CoarseGrid{1+f}(1), obj.FracturesGrid.Grids(f), obj.Coarsening{1+f}(1,:), 1);
                for i=2:obj.maxLevel(1+f)
                    obj.CoarseGrid{1+f}(i).CoarseFactor = obj.Coarsening{1+f}(i,:);
                    obj.CoarseGrid{1+f}(i).BuildCoarseGrid(obj.FracturesGrid.Grids(f));
                    obj.GridMapper.BuildFamily(obj.CoarseGrid{1+f}(i), obj.CoarseGrid{1+f}(i-1), obj.Coarsening{1+f}(1,:), i);
                end

                % Fathers and Verteces
                obj.GridMapper.AssignFathersandVerteces(obj.FracturesGrid.Grids(f), obj.CoarseGrid{1+f}, obj.maxLevel(1+f))
            end
        end
        function CoarseWells(obj, Inj, Prod)
            for i=1:length(Inj)
                % Flag coarse Nodes with wells
                I = Inj(i).Cells;
                for x = 1:obj.maxLevel(1)
                    for j =1:length(I)
                        [r, ~] = find(obj.CoarseGrid{1}(x).GrandChildren == I(j)); % Only in reservoir for now
                        obj.CoarseGrid{1}(x).Wells(r) = 1;
                    end
                end
            end
            for i =1:length(Prod)
                P = Prod(i).Cells;
                for x = 1:obj.maxLevel(1)
                    for j=1:length(P)
                        [r, ~] = find(obj.CoarseGrid{1}(x).GrandChildren == P(j));
                        obj.CoarseGrid{1}(x).Wells(r) = 1;
                    end
                end
            end
        end
        function SelectADMGrid(obj, ProductionSystem)
            % Build ADM Grid for reservoir
            obj.ADMGridSelector.SelectGrid(obj.ReservoirGrid, obj.CoarseGrid{1}, obj.ADMGrid{1}, ProductionSystem, obj.maxLevel);
            obj.ADMStats.N = obj.ADMGrid{1}.N;
            
            for f = 1 : length(obj.maxLevel) - 1
                
            end
        end
        function BuildADMOperators(obj)
            % Build ADM R and P operators
            obj.OperatorsHandler.BuildADMOperators(obj.ReservoirGrid, obj.CoarseGrid{1}, obj.ADMGrid{1});
        end
        function [R, P] = AssembleFullOperators(obj)
            [R, P] = obj.OperatorsHandler.AssembleFullOperators();
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            Formulation.AverageMassOnCoarseBlocks(Status, FluidModel, obj.OperatorsHandler.ADMRest);
        end
    end
end