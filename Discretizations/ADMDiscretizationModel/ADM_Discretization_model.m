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
        ADMGrid_Reservoir
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
            % obj.ADMGrid_Reservoir = adm_grid();
        end
        function AddADMGridSelector(obj, gridselector)
            obj.ADMGridSelector = gridselector;
        end
        function InitializeMapping(obj, ProductionSystem, FluidModel)
            fprintf('Algebraic Dynamic Multilevel (ADM) method run with %d levels and tolerance %s%s = %f\n', ...
                 obj.maxLevel(1), char(916), obj.ADMGridSelector.key, obj.ADMGridSelector.tol);
            % Construct Coarse Grids
            disp(char(2));
            disp('Constructing coarse grids');
            obj.ConstructCoarseGrids(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            obj.FlagPerforatedCoarseCells(ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod);
            
            % Assigning obj.FineGrid and obj.Nf
            if ProductionSystem.FracturesNetwork.Active
                obj.Nf = [obj.ReservoirGrid.N; obj.FracturesGrid.N];
                obj.FineGrid = [obj.ReservoirGrid; obj.FracturesGrid.Grids];

                % Adding the harmonic permeabilities to CrossConnections
                obj.AddHarmonicPermeabilities(ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures);
            else
                obj.FineGrid = obj.ReservoirGrid;
                obj.Nf = obj.ReservoirGrid.N;
            end
            % Global grids (based on global ordering)
            obj.ConstructGlobalGrids();
            
            % Initialise ADM Grid Selector
            obj.ADMGridSelector.Initialise(ProductionSystem, obj.FineGrid, FluidModel.NofPhases);
            
            %% Interpolators
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
                   if obj.CoarseGrid(1,2).Wells{obj.CoarseGrid(1,1).Fathers(i, 1)} == 1
                       obj.ADMGridSelector.NoWellsCoarseCells(i) = 0;
                   end
               end
            else
               for i =1:Nc1
                   if obj.CoarseGrid(1,1).Wells{i} == 1
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
                        for c = 1 : obj.CoarseGrid(1,x).N
                            if any(obj.CoarseGrid(1,x).Children{c,end} == I(j))
                                obj.CoarseGrid(1,x).Wells{c} = 1;
                            end
                        end
                    end
                end
            end
            for i =1:length(Prod)
                P = Prod(i).Cells;
                for x = 1:obj.maxLevel(1)
                    for j=1:length(P)
                        for c = 1 : obj.CoarseGrid(1,x).N
                            if any(obj.CoarseGrid(1,x).Children{c,end} == P(j))
                                obj.CoarseGrid(1,x).Wells{c} = 1;
                            end
                        end
                    end
                end
            end
        end
        function ConstructGlobalGrids(obj)
            %% Create grids based on global ordering
            [n_media, n_levels] = size(obj.CoarseGrid);
            
            %% 1. Matrix of cumulative grid numbers per each media for finesclae and all coarse levels
            Nc_global = zeros(n_media, n_levels+1);
            for m=1:n_media
                Nc_global(m,1) = sum([obj.FineGrid(1:m).N]);
                for L = 1:n_levels
                    Nc_global(m,L+1) = sum([obj.CoarseGrid(1:m,L).N]); 
                end
            end
            
            %% 2. GlobalGrids at fine scale
            obj.GlobalGrids(1).N = sum([obj.FineGrid.N]);
            obj.GlobalGrids(1).Children = cell(obj.GlobalGrids(1).N,1); % at finescale there is no children!
            obj.GlobalGrids(1).Verteces = vertcat(obj.FineGrid(:).Verteces); % This only contains information (zero and one) to check which finescale cell is vertex for higher coarsening levels
            obj.GlobalGrids(1).CoarseFactor = vertcat(obj.FineGrid(:).CoarseFactor);
            obj.GlobalGrids(1).CoarseLevel = 1;
            % 2.1 Fathers
            for x = 1:n_levels - 0
                % 2.1.1 Fathers for the first medium (reservoir)
                Father_Temp = obj.FineGrid(1).Fathers(:,x);
                % 2.1.2 Fathers for the rest of the media (fractures)
                for m = 2:n_media
                    Father_Temp = vertcat(Father_Temp, obj.FineGrid(m).Fathers(:,x) + Nc_global(m-1,x+1));
                end
                obj.GlobalGrids(1).Fathers(:,x) = Father_Temp;
            end
            
            %% 3. GlobalGrids at coarse scales
            for L = 1:n_levels
                obj.GlobalGrids(L+1).N = sum([obj.CoarseGrid(:,L).N]);
                obj.GlobalGrids(L+1).Verteces = vertcat(obj.CoarseGrid(:,L).Verteces); % This only contains information (zero and one) to check which coarsescale cell is vertex for higher coarsening levels
                obj.GlobalGrids(L+1).CoarseFactor = vertcat(obj.CoarseGrid(:,L).CoarseFactor);
                obj.GlobalGrids(L+1).CoarseLevel = L;
                % 3.1 Fathers
                for x = 1:n_levels - L
                    % 3.1.1 Fathers for the first medium (reservoir)
                    Father_Temp = obj.CoarseGrid(1,L).Fathers(:,x);
                    % 3.1.2 Fathers for the rest of the media (fractures)
                    for m = 2:n_media
                        if ~isempty(obj.CoarseGrid(m,L).Fathers)
                            Father_Temp = vertcat(Father_Temp, obj.CoarseGrid(m,L).Fathers(:,x) + Nc_global(m-1,x+2));
                        end
                    end
                    obj.GlobalGrids(L+1).Fathers(:,x) = Father_Temp;
                end
                % 3.2 Children
                for x = 1:L
                    % 3.2.1 Children for the first medium (reservoir)
                    Children_Temp = obj.CoarseGrid(1,L).Children(:,x);
                    % 3.2.2 Children for the rest of the media (fractures)
                    for m = 2:n_media
                        if ~isempty(obj.CoarseGrid(m,L).Children)
                            Children_Temp = vertcat(Children_Temp, cellfun(@(a) a+Nc_global(m-1,x) , obj.CoarseGrid(m,L).Children(:,x),'un',0) );
                        end
                    end
                    obj.GlobalGrids(L+1).Children(:,x) = Children_Temp;
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
        function SelectADMGridCoarse(obj, ProductionSystem, Residual)
            % Build ADM Grid
            obj.ADMGridSelector.SelectGridCoarse(obj.FineGrid, obj.CoarseGrid, obj.ADMGrid, ProductionSystem, Residual, obj.maxLevel);
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
        %% MODIFY PERM:
        function ModifyPerm(obj, ProductionSystem)         % or should I select the ADMgrid?
            ProductionSystem.Reservoir.K = ProductionSystem.Reservoir.K_coarse{1};
            for level = 1:length(obj.CoarseGrid)
                for c =1: obj.CoarseGrid(1, level).N
                    if obj.CoarseGrid(1, level).Active(c) == 1
                        FineCells = cell2mat(obj.CoarseGrid(1, level).GrandChildren(c, :));
                        ProductionSystem.Reservoir.K(FineCells, 1) = ProductionSystem.Reservoir.K_coarse{1 + level}(c, 1);
                        ProductionSystem.Reservoir.K(FineCells, 2) = ProductionSystem.Reservoir.K_coarse{1 + level}(c, 2);
                        ProductionSystem.Reservoir.K(FineCells, 3) = ProductionSystem.Reservoir.K_coarse{1 + level}(c, 3);
                    end
                end
            end
            obj.ReservoirGrid.ComputeRockTransmissibilities(ProductionSystem.Reservoir.K);
        end
        function [Km_Original, Kf_Original] = ModifyPermeabilityContrasts(obj, ProductionSystem)
            % Modifying permeabilities to limit contrast for computation of 
            % coupled basis functions
            
            % Modifying matrix permeability (for now it is only isotropic)
            Km_Original = ProductionSystem.Reservoir.K;
            KmAvg = mean(Km_Original(:));
            if ~isequal( min(Km_Original(:)) , max(Km_Original(:)) )
                K_Log10 = log10(Km_Original(:,1));
                Ratio = (max(K_Log10) - min(K_Log10))/2;
                K_Log10 = ( (K_Log10 - log10(KmAvg)) / Ratio ) + log10(KmAvg);
                ProductionSystem.Reservoir.K(:, 1) = 10.^(K_Log10);
                ProductionSystem.Reservoir.K(:, 2) = 10.^(K_Log10);
                ProductionSystem.Reservoir.K(:, 3) = 10.^(K_Log10);
            end
            KmAvg = mean(ProductionSystem.Reservoir.K(:));
            
            % Modifying fractures permeability
            if ProductionSystem.FracturesNetwork.Active
                dm = mean([obj.ReservoirGrid.dx, obj.ReservoirGrid.dy, obj.ReservoirGrid.dz]);
                Kf_Original = cell(ProductionSystem.FracturesNetwork.NumOfFrac, 1);
                for f=1:ProductionSystem.FracturesNetwork.NumOfFrac
                    Kf_Original{f} = ProductionSystem.FracturesNetwork.Fractures(f).K;
                    KfAvg = mean(Kf_Original{f}(:));
                    if ~isequal( min(Kf_Original{f}(:)) , max(Kf_Original{f}(:)) )
                        ProductionSystem.FracturesNetwork.Fractures(f).K( Kf_Original{f} > KfAvg*10 ) = KfAvg*10;
                        ProductionSystem.FracturesNetwork.Fractures(f).K( Kf_Original{f} < KfAvg/10 ) = KfAvg/10;
                    end
                    KfAvg = mean(ProductionSystem.FracturesNetwork.Fractures(f).K(:));
                    df = obj.FracturesGrid.Grids(f).dz  ;
                    
                    % Reducing the contrast between matrix and fractures
                    Contrast = (KmAvg/dm) ./ (KfAvg/df);
                    Ratio = 1e2;
                    if Contrast > 1
                        multiplier = Contrast / Ratio;
                    else
                        multiplier = Contrast * Ratio;
                    end
                    ProductionSystem.FracturesNetwork.Fractures(f).K = ProductionSystem.FracturesNetwork.Fractures(f).K * multiplier;
                end
            end
        end
        function ResetPermeabilityContrasts(obj, ProductionSystem, Km_Original, Kf_Original)
            % Resetting the permeabilities of matrix and fractures
            ProductionSystem.Reservoir.K = Km_Original;
            if ProductionSystem.FracturesNetwork.Active
                for f=1:ProductionSystem.FracturesNetwork.NumOfFrac
                    ProductionSystem.FracturesNetwork.Fractures(f).K = Kf_Original{f};
                end
                % Adding the harmonic permeabilities to CrossConnections
                obj.AddHarmonicPermeabilities(ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures);
            end
        end
    end
end