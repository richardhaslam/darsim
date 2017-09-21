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
        GridMapper
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
            disp(char(2));
            disp('Constructing coarse grids');
            obj.ConstructCoarseGrids(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            
            if ProductionSystem.FracturesNetwork.Active
                obj.Nf = [obj.ReservoirGrid.N; obj.FracturesGrid.N];
                obj.FineGrid = [obj.ReservoirGrid, obj.FracturesGrid.Grids];
            else
                obj.FineGrid = obj.ReservoirGrid;
                obj.Nf = obj.ReservoirGrid.N;
            end
            % Global grids (based on global ordering)
            obj.ConstructGlobalGrids();
            
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
        function SelectADMGrid(obj, ProductionSystem)
            % Build ADM Grid
            obj.ADMGridSelector.SelectGrid(obj.FineGrid, obj.CoarseGrid, obj.ADMGrid, ProductionSystem, obj.maxLevel);
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