%%  Prolongation builder MS pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 19 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_MSRockTemperature < prolongation_builder
    properties
        C % correction functions
        BFUpdater
        Dimensions
        ADMmap
    end
    methods
        function obj = prolongation_builder_MSRockTemperature(n, cf)
            obj@prolongation_builder(n)
            obj.R = cell(n, 1);
            obj.P = cell(n, 1);
            if cf(1,3) == 1 && cf(1,2) == 1
                obj.Dimensions = 1;
            elseif cf(1,3) == 1
                obj.Dimensions = 2;
            else
                obj.Dimensions = 3;
            end
            obj.ADMmap = adm_map(prod(cf, 2));
        end
        function BuildStaticOperators(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections, maxLevel, CoarseGrid)
            % Initialise
            obj.R = cell(maxLevel(1), 1);
            obj.P = cell(maxLevel(1), 1);
            % Build Rock Temperaature system
            obj.BFUpdater.ConstructRockTemperatureSystem(ProductionSystem, FineGrid, CrossConnections);
            %Build static restriction operator (FV)
            disp('Building Rock Temperature Restriction 1');
            start1 = tic;
            obj.R{1} = obj.MsRestriction(FineGrid(1), CoarseGrid(1,1));
            % Build Prolongation operator
            disp('Building Rock Temperature Prolongation 1');
            [obj.P{1}, obj.C{1}] = obj.BFUpdater.MsProlongation(FineGrid(1), CoarseGrid(1,1), obj.Dimensions);
            % Build tpfa coarse system of level 1 (with MsFE)
            obj.BFUpdater.UpdatePressureMatrix(obj.P{1}, CoarseGrid(1, 1));
            for x = 2:maxLevel(1)
                % Build static restriction operator (FV)
                disp(['Building Rock Temperature Restriction ', num2str(x)]);
                obj.R{x} = obj.MsRestriction(CoarseGrid(1, x-1), CoarseGrid(1, x));
                % Build Prolongation operator
                disp(['Building Rock Temperature Prolongation ', num2str(x)]);
                [obj.P{x}, obj.C{x}] = obj.BFUpdater.MsProlongation(CoarseGrid(1, x-1), CoarseGrid(1, x), obj.Dimensions);
                %Build tpfa coarse system of level x (with MsFE)
                obj.BFUpdater.UpdatePressureMatrix(obj.P{x}, CoarseGrid(1, x));
            end
            StaticOperators = toc(start1);
            disp(['Temperature static operators built in: ', num2str(StaticOperators), ' s']);
        end
        function ADMProlTr = ADMProlongation(obj, ADMGrid_original, GlobalGrids, ADMRest)
            start1 = tic;
            %Copy ADM Grid
            ADMGrid = ADMGrid_original.copy();
            
            % Modify ADM_Grid to remove fratures
            ADMGrid_Reservoir = obj.ExtractReservoirADMGrid(ADMGrid);
%             ADMGrid_Reservoir = ADMGrid;
                     
            % Pressure prolongation
            ADMProlTr = 1;
            
            % Loop over the levels
            for level = ADMGrid_Reservoir.MaxLevel:-1:2
                ProlTr = obj.LevelProlongation(ADMGrid_Reservoir, GlobalGrids(level), level);
                % Multiply by previous objects
                ADMProlTr = ProlTr * ADMProlTr;
            end
            
            % Last prolongation is different coz I use fine-scale ordering
            ProlTr = obj.LastProlongation(ADMGrid_Reservoir, GlobalGrids(1));
            
            % Multiply by previous objects
            ADMProlTr = ProlTr * ADMProlTr;
            
            prolongation = toc(start1);
            disp(['Rock Temperature Prolongation built in: ', num2str(prolongation), ' s']);
        end
        function ProlTr = LevelProlongation(obj, ADMGrid, FineGrid, level)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For a given level the prolongation operator looks like this
            %       Nf     Nc
            %    ----       ----
            %    |      |      |
            % Nf |  I   |  0   |
            %    |______|______|
            %    |      |      |
            % Nx |      |      |
            %    |      |      |
            %    ----       ----
            
            
            % Update map for the next level
            obj.ADMmap.Update(ADMGrid, FineGrid, level);
            
            % 1. Build the object
            ProlTr = sparse(sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nx), sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nc));
            
            % 2. Fill in top left
            ProlTr(1:sum(obj.ADMmap.Nf), 1:sum(obj.ADMmap.Nf)) = speye(sum(obj.ADMmap.Nf));
            
            % 3. Fill in Bottom left
            ProlTr(sum(obj.ADMmap.Nf) + 1 : end,  obj.ADMmap.Verteces) = obj.P{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexVerteces);
            
            % 4. Fill in Bottom right
            ProlTr(sum(obj.ADMmap.Nf) + 1 :end, sum(obj.ADMmap.Nf) + 1 : end) = obj.P{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexNc);
        end
        function ProlTr = LastProlongation(obj, ADMGrid, FineGrid)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For the first level the prolongation operator looks like this
            %       Nf     Nl1
            %    ----       ----
            %    |      |      |
            %    |      |   0  |
            % Nl0|      |______|
            %    |      |      |
            %    |      |      |
            %    |      |      |
            %    ----       ----

            % Update map
            obj.ADMmap.Update(ADMGrid, FineGrid, 1);
       
            % 1. Build object
            ProlTr = sparse(sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nx), sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nc));
            
            % 2. Fill in FS verteces of level 1
            ProlTr(:,  obj.ADMmap.Verteces) = obj.P{1}(:, obj.ADMmap.OriginalIndexVerteces);
            
            % 3. Fill in coarse-scale nodes
            ProlTr(:, sum(obj.ADMmap.Nf) + 1 : end) = obj.P{1}(:, obj.ADMmap.OriginalIndexNc);
            
            % 4. Fill in fine-scale nodes first
            rows = obj.ADMmap.OriginalIndexNf';
            columns = 1:sum(obj.ADMmap.Nf);
            ProlTr(rows, :) = 0; % if it s fine-scale already I get rid of useless fillings
            ProlTr(sub2ind(size(ProlTr), rows, columns)) = 1;
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            % for now I do not update pressure basis functions
        end
        function AverageMassOnCoarseBlocks(obj, Formulation, ProductionSystem, FineGrid, FluidModel, ADMRest)
            % virtual call
        end
        function ADMGrid_Reservoir = ExtractReservoirADMGrid(obj,ADMGrid)
            % Remove the fractures from ADMGrid to be used for RockTemperature Prolongation
            ADMGrid_Reservoir = adm_grid();
            ADMGrid_Reservoir.N = ADMGrid.N(1,:);
            ADMGrid_Reservoir.Ntot = sum(ADMGrid_Reservoir.N);
            ADMGrid_Reservoir.MaxLevel = ADMGrid.MaxLevel;
            ADMGrid_Reservoir.level = zeros(ADMGrid_Reservoir.Ntot,1);
            ADMGrid_Reservoir.CellIndex = zeros(ADMGrid_Reservoir.Ntot,1);
            ADMGrid_Reservoir.Children = cell(ADMGrid_Reservoir.Ntot,1);
            ADMGrid_Reservoir.GrandChildren = cell(ADMGrid_Reservoir.Ntot,1);
            for level = 0 : length(ADMGrid_Reservoir.N)-1
                IndexResStrart  = sum(ADMGrid_Reservoir.N(1:level  )) + 1;
                IndexResEnd     = sum(ADMGrid_Reservoir.N(1:level+1));
                IndexFullStrart = sum(sum(ADMGrid.N(:,1:level))) + 1 ;
                IndexFullEnd    = sum(sum(ADMGrid.N(:,1:level))) + ADMGrid.N(1,level+1);
                ADMGrid_Reservoir.level        ( IndexResStrart : IndexResEnd , 1 ) = ADMGrid.level        ( IndexFullStrart : IndexFullEnd , 1 );
                ADMGrid_Reservoir.CellIndex    ( IndexResStrart : IndexResEnd , 1 ) = ADMGrid.CellIndex    ( IndexFullStrart : IndexFullEnd , 1 );
                ADMGrid_Reservoir.Fathers      ( IndexResStrart : IndexResEnd , : ) = ADMGrid.Fathers      ( IndexFullStrart : IndexFullEnd , : );
                ADMGrid_Reservoir.Children     ( IndexResStrart : IndexResEnd , 1 ) = ADMGrid.Children     ( IndexFullStrart : IndexFullEnd , 1 );
                ADMGrid_Reservoir.GrandChildren( IndexResStrart : IndexResEnd , 1 ) = ADMGrid.GrandChildren( IndexFullStrart : IndexFullEnd , 1 );
                ADMGrid_Reservoir.Verteces     ( IndexResStrart : IndexResEnd , : ) = ADMGrid.Verteces     ( IndexFullStrart : IndexFullEnd , : );
                ADMGrid_Reservoir.CoarseFactor ( IndexResStrart : IndexResEnd , : ) = ADMGrid.CoarseFactor ( IndexFullStrart : IndexFullEnd , : );
            end    
        end
    end
end