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
            
            % Modifying permeabilities to limit contrast for computation of 
            % coupled basis functions
            
            % Modifying matrix permeability (for now isotropic)
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
                        ProductionSystem.FracturesNetwork.Fractures(f).K( Kf_Original > KfAvg*10 ) = KfAvg*10;
                        ProductionSystem.FracturesNetwork.Fractures(f).K( Kf_Original < KfAvg/10 ) = KfAvg/10;
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

            % Assigning obj.FineGrid and obj.Nf
            if ProductionSystem.FracturesNetwork.Active
                obj.Nf = [obj.ReservoirGrid.N; obj.FracturesGrid.N];
                obj.FineGrid = [obj.ReservoirGrid, obj.FracturesGrid.Grids];

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
            
            % Resetting the permeabilities of matrix and fractures
            ProductionSystem.Reservoir.K = Km_Original;
            if ProductionSystem.FracturesNetwork.Active
                for f=1:ProductionSystem.FracturesNetwork.NumOfFrac
                    ProductionSystem.FracturesNetwork.Fractures(f).K = Kf_Original{f};
                end
                % Adding the harmonic permeabilities to CrossConnections
                obj.AddHarmonicPermeabilities(ProductionSystem.Reservoir, ProductionSystem.FracturesNetwork.Fractures);
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
                            % [r, ~] = find([obj.CoarseGrid(1,x).GrandChildren{c,:}] == I(j)); % Only in reservoir for now
                            if any(obj.CoarseGrid(1,x).GrandChildren{c,:} == I(j))
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
                            % [r, ~] = find([obj.CoarseGrid(1,x).GrandChildren{c,:}] == P(j));
                            if any(obj.CoarseGrid(1,x).GrandChildren{c,:} == P(j))
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
            if any(strcmp(ProductionSystem.Reservoir.State.Properties.keys,'Tr'))
                obj.ExtractReservoirADMGrid;
                obj.OperatorsHandler.BuildReservoirADMRestriction(obj.ADMGrid_Reservoir, obj.ReservoirGrid)
            end
            obj.OperatorsHandler.BuildADMOperators(obj.GlobalGrids, obj.ADMGrid);
        end
        function SelectADMGridCoarse(obj, ProductionSystem, Residual)
            % Build ADM Grid
            obj.ADMGridSelector.SelectGridCoarse(obj.FineGrid, obj.CoarseGrid, obj.ADMGrid, ProductionSystem, Residual, obj.maxLevel);
            obj.ADMStats.N = obj.ADMGrid.N(1,:);
            
            % Update prolongation builders
            obj.OperatorsHandler.UpdateProlongationOperators(obj.FineGrid, obj.CoarseGrid, ProductionSystem);
            
            % Build ADM R and P operators
            if any(strcmp(ProductionSystem.Reservoir.State.Properties.keys,'Tr'))
                obj.ExtractReservoirADMGrid;
                obj.OperatorsHandler.BuildReservoirADMRestriction(obj.ADMGrid_Reservoir, obj.ReservoirGrid)
            end
            obj.OperatorsHandler.BuildADMOperators(obj.GlobalGrids, obj.ADMGrid);
        end
        function ExtractReservoirADMGrid(obj)
            % Remove the fractures from ADMGrid to be used for RockTemperature Prolongation
            obj.ADMGrid_Reservoir = adm_grid();
            obj.ADMGrid_Reservoir.N = obj.ADMGrid.N(1,:);
            obj.ADMGrid_Reservoir.Ntot = sum(obj.ADMGrid_Reservoir.N);
            obj.ADMGrid_Reservoir.MaxLevel = obj.ADMGrid.MaxLevel;
            obj.ADMGrid_Reservoir.level = zeros(obj.ADMGrid_Reservoir.Ntot,1);
            obj.ADMGrid_Reservoir.CellIndex = zeros(obj.ADMGrid_Reservoir.Ntot,1);
            obj.ADMGrid_Reservoir.Children = cell(obj.ADMGrid_Reservoir.Ntot,1);
            obj.ADMGrid_Reservoir.GrandChildren = cell(obj.ADMGrid_Reservoir.Ntot,1);
            for level = 0 : length(obj.ADMGrid_Reservoir.N)-1
                IndexResStrart  = sum(obj.ADMGrid_Reservoir.N(1:level  )) + 1;
                IndexResEnd     = sum(obj.ADMGrid_Reservoir.N(1:level+1));
                IndexFullStrart = sum(sum(obj.ADMGrid.N(:,1:level))) + 1 ;
                IndexFullEnd    = sum(sum(obj.ADMGrid.N(:,1:level))) + obj.ADMGrid.N(1,level+1);
                obj.ADMGrid_Reservoir.level        ( IndexResStrart : IndexResEnd , 1 ) = obj.ADMGrid.level        ( IndexFullStrart : IndexFullEnd , 1 );
                obj.ADMGrid_Reservoir.CellIndex    ( IndexResStrart : IndexResEnd , 1 ) = obj.ADMGrid.CellIndex    ( IndexFullStrart : IndexFullEnd , 1 );
                obj.ADMGrid_Reservoir.Fathers      ( IndexResStrart : IndexResEnd , : ) = obj.ADMGrid.Fathers      ( IndexFullStrart : IndexFullEnd , : );
                obj.ADMGrid_Reservoir.Children     ( IndexResStrart : IndexResEnd , 1 ) = obj.ADMGrid.Children     ( IndexFullStrart : IndexFullEnd , 1 );
                obj.ADMGrid_Reservoir.GrandChildren( IndexResStrart : IndexResEnd , 1 ) = obj.ADMGrid.GrandChildren( IndexFullStrart : IndexFullEnd , 1 );
                obj.ADMGrid_Reservoir.Verteces     ( IndexResStrart : IndexResEnd , : ) = obj.ADMGrid.Verteces     ( IndexFullStrart : IndexFullEnd , : );
                obj.ADMGrid_Reservoir.CoarseFactor ( IndexResStrart : IndexResEnd , : ) = obj.ADMGrid.CoarseFactor ( IndexFullStrart : IndexFullEnd , : );

            end    
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
    end
end