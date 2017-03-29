%  ADM discretization model
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
        function obj = ADM_Discretization_model(nx, ny, nz, maxlevel, coarsening)
            obj@Discretization_model(nx, ny, nz);
            obj.CoarseGrid = coarse_grid();
            obj.ADMGrid = adm_grid();
            obj.Coarsening = coarsening;
            obj.maxLevel = maxlevel;
            obj.GridMapper = grid_mapper();
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
            obj.ConstructCoarseGrids(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            %% Pressure interpolators
            disp('Static operators - start computation');
            obj.OperatorsHandler.BuildStaticOperators(obj.CoarseGrid, obj.ReservoirGrid, obj.maxLevel,...
                ProductionSystem.Reservoir.K, ProductionSystem.Reservoir.State.S, FluidModel);
            disp('Static operators - end')
            disp(char(2));
        end
        function ConstructCoarseGrids(obj, Inj, Prod)
            % Construct all coarse grids
            obj.CoarseGrid(1).CoarseFactor = obj.Coarsening(1,:);
            obj.CoarseGrid(1).BuildCoarseGrid(obj.ReservoirGrid);
            obj.GridMapper.BuildFamily(obj.CoarseGrid(1), obj.ReservoirGrid, obj.Coarsening(1,:), 1);
            for i=2:obj.maxLevel
                obj.CoarseGrid(i).CoarseFactor = obj.Coarsening(i,:);
                obj.CoarseGrid(i).BuildCoarseGrid(obj.ReservoirGrid);
                obj.GridMapper.BuildFamily(obj.CoarseGrid(i), obj.CoarseGrid(i-1), obj.Coarsening(1,:), i);
            end
            
            % Fathers and Verteces
            obj.GridMapper.AssignFathersandVerteces(obj.ReservoirGrid, obj.CoarseGrid, obj.maxLevel)

            % Flag coarse blocks with wells
            obj.CoarseWells(Inj, Prod);
            obj.ADMGridSelector.NoWellsCoarseCells = ones(obj.CoarseGrid(1).N, 1);
            Nc1 = obj.CoarseGrid(1).N;
            if obj.maxLevel > 1
               for i = 1:Nc1
                   if obj.CoarseGrid(2).Wells(obj.CoarseGrid(1).Fathers(i, 2)) == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            else
               for i =1:Nc1
                   if obj.CoarseGrid(1).Wells(i) == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            end
        end
        function CoarseWells(obj, Inj, Prod)
            for i=1:length(Inj)
                % Flag coarse Nodes with wells
                I = Inj(i).Cells;
                for x = 1:obj.maxLevel
                    for j =1:length(I)
                        [r, ~] = find(obj.CoarseGrid(x).GrandChildren == I(j));
                        obj.CoarseGrid(x).Wells(r) = 1;
                    end
                end
            end
            for i =1:length(Prod)
                P = Prod(i).Cells;
                for x = 1:obj.maxLevel
                    for j=1:length(P)
                        [r, ~] = find(obj.CoarseGrid(x).GrandChildren == P(j));
                        obj.CoarseGrid(x).Wells(r) = 1;
                    end
                end
            end
        end
        function SelectADMGrid(obj, ProductionSystem)
            % Build ADM Grid
            obj.ADMGridSelector.SelectGrid(obj.ReservoirGrid, obj.CoarseGrid, obj.ADMGrid, ProductionSystem, obj.maxLevel);
            obj.ADMStats.N = obj.ADMGrid.N;
        end
        function BuildADMOperators(obj)
            % Build ADM R and P operators
            obj.OperatorsHandler.BuildADMOperators(obj.ReservoirGrid, obj.CoarseGrid, obj.ADMGrid);
        end
        function AverageMassOnCoarseBlocks(obj, Status, FluidModel, Formulation)
            Formulation.AverageMassOnCoarseBlocks(Status, FluidModel, obj.OperatorsHandler.ADMRest, obj.OperatorsHandler.ADMProlAv);
        end
    end
end