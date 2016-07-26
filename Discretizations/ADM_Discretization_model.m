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
        ADMSettings
        CoarseGrid
        ADMGrid
    end
    methods
        function obj = ADM_Discretization_model(nx, ny, nz, settings)
            obj@Discretization_model(nx, ny, nz);
            obj.ADMSettings = settings;
        end
        function Initialize(obj, ProductionSystem)
            disp('Algebraic Dynamic Multilevel (ADM) method run')
            obj.ReservoirGrid.Initialize(ProductionSystem.Reservoir);
            
            % Construct Coarse Grids
            obj.ConstructCoarseGrids();
            
            %% Pressure interpolators
            disp('Pressure interpolator - start computation');
            obj.StaticPressureInterpolators(ProductionSystem.Reservoir.K);
            disp('Pressure interpolator - end')
            disp(char(2));
        end
        function ConstructCoarseGrids(obj)
            obj.ReservoirGrid.AddCoordinates();
            
            %Build Coarse Grids
            CoarseGrid = struct('CoarseFactor', {}, 'Nx', {}, 'Ny', {}, ...
                'I', {}, 'J', {},'Father', {}, 'Active', {}, 'Wells', {}, 'Neighbours', {}, 'Centers', {});
            for i=1:obj.ADMSettings.maxLevel
                obj.CoarseGrid(i) = coarse_grid();
            end
            [Grid, CoarseGrid] = AssignFathers(Grid, CoarseGrid, ADMSettings.maxLevel);
            
            %Flag coarse blocks with wells
            [CoarseGrid] = CoarseWells(Grid, CoarseGrid, maxLevel, Inj, Prod);
            NoWellsCoarseCells = ones(CoarseGrid(1).Nx*CoarseGrid(1).Ny,1);
            Nc1 = CoarseGrid(1).Nx * CoarseGrid(1).Ny;
            if maxLevel > 1
                for i = 1:Nc1
                    if CoarseGrid(2).Wells(CoarseGrid(1).Father(i,2)) == 1
                        NoWellsCoarseCells(i) = 0;
                    end
                end
            else
                for i =1:Nc1
                    if CoarseGrid(1).Wells(i) == 1
                        NoWellsCoarseCells(i) = 0;
                    end
                end
            end
        end
    end
end