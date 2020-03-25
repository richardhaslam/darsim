%  MS discretization model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Multiscale_Discretization_model < Discretization_model
    properties
        Nf
        Nc
        Vertex_On_Corner
        maxLevel
        CoarseningSwitch % This is a switch for Coarsening options (mostly valid for fractures).
        % "0" means coarsened upto the max coarsening level, afterwards stays at that resolution (higher computational demand)
        % "1" means coarsened upto the max coarsening level, afterwards no coare node (lower computational demand)
        Coarsening
        CoarseGrid
        OperatorsHandler
        FineGrid
        GridMapper
    end
    methods
        function obj = Multiscale_Discretization_model(maxlevel, coarsening)
            n = length(maxlevel);
            obj.CoarseGrid = coarse_grid.empty;
            obj.Coarsening = coarsening;
            obj.maxLevel   = maxlevel;
            obj.Nf = zeros(n, 1);
            obj.Nc = zeros(n, maxlevel(1));
            obj.GridMapper = grid_mapper();
        end
        function AddOperatorsHandler(obj, operatorshandler)
            obj.OperatorsHandler = operatorshandler;
        end
        function InitializeMapping(obj, ProductionSystem, FluidModel)
            disp([num2str(obj.maxLevel(1)), ' levels multiscale run']);
            % Construct Coarse Grids
            disp(newline);
            disp('Constructing coarse grids');
            obj.ConstructCoarseGrids(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            if ProductionSystem.FracturesNetwork.Active
                obj.Nf = [obj.ReservoirGrid.N; obj.FracturesGrid.N];
                obj.FineGrid = [obj.ReservoirGrid, obj.FracturesGrid.Grids];
            else
                obj.FineGrid = obj.ReservoirGrid;
                obj.Nf = obj.ReservoirGrid.N;
            end

            %% Pressure interpolators
            disp('Multiscale basis functions - start computation');
            start = tic;
            obj.OperatorsHandler.BuildMMsOperators(ProductionSystem, FluidModel, obj.FineGrid, obj.CrossConnections, ...
                obj.maxLevel, obj.CoarseGrid);
            disp('Multiscale basis functions - end')
            timer = toc(start);
            disp(['Multiscale basis functions computation took ', num2str(timer)])
            disp(newline);
        end
        function ConstructCoarseGrids(obj, Inj, Prod)
            if mod( obj.ReservoirGrid.Nx , obj.Coarsening(1,1,1) ) == 0 && ...
               mod( obj.ReservoirGrid.Ny , obj.Coarsening(1,2,1) ) == 0 && ...
               mod( obj.ReservoirGrid.Nz , obj.Coarsening(1,3,1) ) == 0
                obj.Vertex_On_Corner = 0;
                fprintf('The verteces are not on the corners.\n');
            elseif mod( max((obj.ReservoirGrid.Nx-1),1) , obj.Coarsening(1,1,1) ) == 0 && ...
                   mod( max((obj.ReservoirGrid.Ny-1),1) , obj.Coarsening(1,2,1) ) == 0 && ...
                   mod( max((obj.ReservoirGrid.Nz-1),1) , obj.Coarsening(1,3,1) ) == 0
                obj.Vertex_On_Corner = 1;
                fprintf('The verteces are on the corners.\n');
            else
                error('The number of reservoir grid cells do not comply with coarsening ratios. Please check the input files.');
            end
            %% 1. Reservoir
            % Construct all coarse grids for reservoir
            obj.CoarseGrid(1,1) = coarse_grid();
            obj.CoarseGrid(1,1).CoarseFactor = obj.Coarsening(1,:,1);
            obj.CoarseGrid(1,1).Vertex_On_Corner = obj.Vertex_On_Corner;
            obj.CoarseGrid(1,1).BuildCoarseGrid(obj.ReservoirGrid, obj.Coarsening(1,:,1));
            obj.GridMapper.BuildFamily(obj.CoarseGrid(1,1), obj.ReservoirGrid, obj.Coarsening(1,:,1), 1);
            obj.CoarseGrid(1,1).AddWells(Inj, Prod);
            obj.Nc(1, 1) = obj.CoarseGrid(1,1).N;
            for L=2:obj.maxLevel(1)
                obj.CoarseGrid(1,L) = coarse_grid();
                obj.CoarseGrid(1,L).CoarseFactor = obj.Coarsening(1,:,L);
                obj.CoarseGrid(1,L).Vertex_On_Corner = obj.Vertex_On_Corner;
                obj.CoarseGrid(1,L).BuildCoarseGrid(obj.CoarseGrid(1,L-1),obj.Coarsening(1,:,L)./obj.Coarsening(1,:,L-1));
                %obj.CoarseGrid(1,L).BuildCoarseGrid(obj.ReservoirGrid);
                obj.GridMapper.BuildFamily(obj.CoarseGrid(1,L), obj.CoarseGrid(1,L-1), obj.Coarsening(1,:,L)./obj.Coarsening(1,:,L-1), L);
                obj.CoarseGrid(1,L).AddWells(Inj, Prod);
                obj.Nc(1, L) = obj.CoarseGrid(1,L).N;
            end
            
            % Fathers and Verteces for reservoir
            obj.GridMapper.AssignFathersandVerteces(obj.ReservoirGrid, obj.CoarseGrid(1,:), obj.maxLevel(1))
            
            %% 2. Fractures
            % Construct all coarse grids for fractures
            for f = 1 : size(obj.Coarsening,1) - 1
                % minMaxLevel = min( obj.maxLevel(1) , obj.maxLevel(1+f) );
                obj.CoarseGrid(1+f,1) = coarse_grid();
                obj.CoarseGrid(1+f,1).CoarseFactor = obj.Coarsening(1+f,:,1);
                obj.CoarseGrid(1+f,1).Vertex_On_Corner = obj.Vertex_On_Corner;
                obj.CoarseGrid(1+f,1).BuildCoarseGrid(obj.FracturesGrid.Grids(f),obj.Coarsening(1+f,:,1));
                obj.GridMapper.BuildFamily(obj.CoarseGrid(1+f,1), obj.FracturesGrid.Grids(f), obj.Coarsening(1+f,:,1), 1);
                obj.Nc(f+1,1) = obj.CoarseGrid(1+f,1).N;
                for L = 2 : obj.maxLevel(1)
                    obj.CoarseGrid(1+f,L) = coarse_grid();
                    obj.CoarseGrid(1+f,L).CoarseFactor = obj.Coarsening(1+f,:,L);
                    obj.CoarseGrid(1+f,L).Vertex_On_Corner = obj.Vertex_On_Corner;
                    obj.CoarseGrid(1+f,L).BuildCoarseGrid(obj.CoarseGrid(1+f,L-1),obj.Coarsening(1+f,:,L)./obj.Coarsening(1+f,:,L-1));
                    %obj.CoarseGrid(1+f,L).BuildCoarseGrid(obj.FracturesGrid.Grids(f));
                    obj.GridMapper.BuildFamily(obj.CoarseGrid(1+f,L), obj.CoarseGrid(1+f,L-1), obj.Coarsening(1+f,:,L)./obj.Coarsening(1+f,:,L-1), L);
                    obj.Nc(f+1,L) = obj.CoarseGrid(1+f,L).N;
                end
                obj.GridMapper.AssignFathersandVerteces(obj.FracturesGrid.Grids(f), obj.CoarseGrid(1+f,:), obj.maxLevel(1));

                if obj.CoarseningSwitch(1+f) == 1
                    % This fracture is coarsened upto the last coarsening level,
                    % afterwards no coare node (lowest computational demand).
                    for L = obj.maxLevel(1+f)+1 : obj.maxLevel(1)
                        obj.CoarseGrid(1+f,L) = coarse_grid();
                        obj.CoarseGrid(1+f,L).CoarseFactor = obj.Coarsening(1+f,:,L);
                        obj.CoarseGrid(1+f,L).Nx = 0;
                        obj.CoarseGrid(1+f,L).Ny = 0;
                        obj.CoarseGrid(1+f,L).Nz = 0;
                        obj.CoarseGrid(1+f,L).N = 0;
                    end
                end
            end
            
			fprintf('Coarsening ratio in reservoir: %d x %d x %d\n' , obj.Coarsening(1,1,1), obj.Coarsening(1,2,1), obj.Coarsening(1,3,1) );
            for L = 1 : obj.maxLevel(1)
                fprintf('Number of reservoir coarse nodes at level %d: %d\n' , L, obj.Nc(1,L) );
                if (size(obj.Coarsening,1) - 1) > 0
                    fprintf('Number of fractures coarse nodes at level %d: %d\n' , L, sum(obj.Nc(2:end,L)));
                end
            end
        end
        function AddWellsToInitialPressure(obj, ProductionSystem, FluidModel)
            % Improving the first pressure guess for multilevel/multiscale method
            deltaP_w = obj.OperatorsHandler.ProlongationBuilders(1).StaticMultilevelPressureGuess(ProductionSystem, FluidModel, obj.FineGrid, obj.CoarseGrid(:, end), obj.CrossConnections);

            %% 1. Updating Reservoir Pressure with Well Correction
            Pm = ProductionSystem.Reservoir.State.Properties(['P_', num2str(FluidModel.NofPhases)]);
            Pm.update(deltaP_w(1:obj.ReservoirGrid.N));
            
            %% 2. Update fractures pressure and densities
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    % Update Pressure
                    index1 = obj.ReservoirGrid.N + sum(obj.FracturesGrid.N(1:f-1)) + 1;
                    index2 = index1 - 1 + obj.FracturesGrid.N(f);
                    Pf = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(['P_', num2str(FluidModel.NofPhases)]);
                    Pf.update(deltaP_w(index1:index2));
                end
            end
            
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end
        function SelectADMGrid(obj, ProductionSystem)
            % virtual call
        end
    end
end