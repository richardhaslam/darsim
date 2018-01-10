%  ADM discretization model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 7 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef ADM_Discretization_model < Multiscale_Discretization_model
    properties
        ADMGrid
        ADMGridSelector
        ADMStats
        GlobalGrids
    end
    methods
        function obj = ADM_Discretization_model(maxlevel, coarsening)
            % build a Multiscale discretization first
            obj@Multiscale_Discretization_model(maxlevel, coarsening)
            obj.GlobalGrids = grid_darsim.empty;
            obj.ADMGrid    = adm_grid();
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
            obj.FlagPerforatedCoarseCells(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            if ProductionSystem.FracturesNetwork.Active
                obj.Nf = [obj.ReservoirGrid.N; obj.FracturesGrid.N];
                obj.FineGrid = [obj.ReservoirGrid, obj.FracturesGrid.Grids];
            else
                obj.FineGrid = obj.ReservoirGrid;
                obj.Nf = obj.ReservoirGrid.N;
            end
            % Global grids (based on global ordering)
            obj.ConstructGlobalGrids();
            
            % Initialise ADM Grid Selector
            obj.ADMGridSelector.Initialise(ProductionSystem, obj.FineGrid, FluidModel.NofPhases);
            
            %% Pressure interpolators
            disp('Static operators - start computation');
            start = tic;
            for i=1:length(obj.OperatorsHandler.ProlongationBuilders)
                obj.OperatorsHandler.ProlongationBuilders(i).BuildStaticOperators(ProductionSystem, FluidModel, obj.FineGrid, obj.CrossConnections, ...
                    obj.maxLevel, obj.CoarseGrid);
            end
            disp('Static operators - end')
            timer = toc(start);
            disp(['Static operators construction took ', num2str(timer)])
            disp(char(2));
        end
        function FlagPerforatedCoarseCells(obj, Inj, Prod)
            % Flag coarse blocks with wells for reservoir
            obj.CoarseWells(Inj, Prod);
            obj.ADMGridSelector.NoWellsCoarseCells = ones(obj.CoarseGrid(1,1).N, 1);
            Nc1 = obj.CoarseGrid(1,1).N;
            if obj.maxLevel(1) > 1
               for i = 1:Nc1
                   if obj.CoarseGrid(1,2).Wells{obj.CoarseGrid(1,1).Fathers(i, 2)} == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            else
               for i =1:Nc1
                   if obj.CoarseGrid(1,1).Wells(i) == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            end
        end
        function CoarseWells(obj, Inj, Prod)
            for i=1:length(Inj)
                % Flag coarse Nodes with wells
                I = Inj(i).Cells;
                for x = 1:obj.maxLevel(1)
                    for j =1:length(I)
                        [r, ~] = find(obj.CoarseGrid(1,x).GrandChildren == I(j)); % Only in reservoir for now
                        obj.CoarseGrid(1,x).Wells{r} = 1;
                    end
                end
            end
            for i =1:length(Prod)
                P = Prod(i).Cells;
                for x = 1:obj.maxLevel(1)
                    for j=1:length(P)
                        [r, ~] = find(obj.CoarseGrid(1,x).GrandChildren == P(j));
                        obj.CoarseGrid(1,x).Wells{r} = 1;
                    end
                end
            end
        end
        function ConstructGlobalGrids(obj)
            %% Create grids based on global ordering
            [n_media, n_levels] = size(obj.CoarseGrid);
           
            % 1. Initialise global grids
            Nc_global = zeros(n_media, n_levels);
            obj.GlobalGrids(1).N = sum([obj.FineGrid.N]);
            obj.GlobalGrids(1).Initialise(n_levels);
            for m=1:n_media
                Nc_global(m, 1) = sum([obj.FineGrid(1:m-1).N]);
            end
            for x = 1:n_levels
                obj.GlobalGrids(x+1).N = sum([obj.CoarseGrid(:, x).N]);
                obj.GlobalGrids(x+1).Initialise(n_levels);
                for m=1:n_media
                    Nc_global(m, x+1) = sum([obj.CoarseGrid(1:m-1, x).N]);    
                end
            end
            
            for m=1:n_media
                % Global Fine Grid
                obj.GlobalGrids(1).CopyGridEntries(obj.FineGrid(m), Nc_global(m, :), 1);
                % Global Coarse Grid
                for i=1:n_levels
                    obj.GlobalGrids(i+1).CopyGridEntries(obj.CoarseGrid(m, i), Nc_global(m, :), i+1);
                end
            end
        end
        function SelectADMGrid(obj, ProductionSystem, Residual)
            % Build ADM Grid
            obj.ADMGridSelector.SelectGrid(obj.FineGrid, obj.CoarseGrid, obj.ADMGrid, ProductionSystem, Residual, obj.maxLevel);
            obj.ADMStats.N = obj.ADMGrid.N(1,:);
            
            % Update prolongation builders
            obj.OperatorsHandler.UpdateProlongationOperators(obj.FineGrid, obj.CoarseGrid, ProductionSystem);
            
            % Build ADM R and P operators
            obj.OperatorsHandler.BuildADMOperators(obj.GlobalGrids, obj.ADMGrid);
        end
        function [R, P] = AssembleFullOperators(obj)
            [R, P] = obj.OperatorsHandler.AssembleFullOperators();
        end
        function AverageMassOnCoarseBlocks(obj, ProductionSystem, FluidModel, Formulation)
            obj.OperatorsHandler.ProlongationBuilders(2).AverageMassOnCoarseBlocks(Formulation, ProductionSystem, obj.FineGrid, FluidModel, obj.OperatorsHandler.ADMRest);  
        end
    end
end