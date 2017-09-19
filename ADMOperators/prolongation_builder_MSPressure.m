%  Prolongation builder MS pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 August 2017
%Last modified: 19 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef prolongation_builder_MSPressure < prolongation_builder
    properties
        BFUpdater
        Dimensions
        ADMmap
    end
    methods
        function obj = prolongation_builder_MSPressure(n, cf)
            obj@prolongation_builder(n)
            obj.R = cell(1, n);
            obj.P = cell(1, n);
            if cf(3) == 1 && cf(2) == 1
                obj.Dimensions = 1;
            elseif cf(3) == 1
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
            % Build Pressure system
            obj.BFUpdater.ConstructPressureSystem(ProductionSystem, FluidModel, FineGrid, CrossConnections);
            
            %Build static restriction operator (FV)
            disp('Building Restriction 1');
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(:,1));
            % Build Prolongation operator
            disp('Building Prolongation 1');
            obj.P{1} = obj.BFUpdater.MsProlongation(FineGrid, CoarseGrid(:,1), obj.Dimensions);
            % Build tpfa coarse system of level 1 (with MsFE)
            obj.BFUpdater.UpdatePressureMatrix(obj.P{1}, CoarseGrid(:, 1));
            for x = 2:maxLevel(1)
                % Build static restriction operator (FV)
                disp(['Building Restriction ', num2str(x)]);
                obj.R{x} = obj.MsRestriction(CoarseGrid(:, x-1), CoarseGrid(:, x));
                % Build Prolongation operator
                disp(['Building Prolongation ', num2str(x)]);
                obj.P{x} = obj.BFUpdater.MsProlongation(CoarseGrid(:, x-1), CoarseGrid(:, x), obj.Dimensions);
                %Build tpfa coarse system of level x (with MsFE)
                obj.BFUpdater.UpdatePressureMatrix(obj.P{x}, CoarseGrid(:, x));
            end
        end
        function ADMProlp = ADMProlongation(obj, ADMGrid, GlobalGrids, ADMRest)
            % Pressure prolongation
            ADMProlp = 1;
            
            % Loop over the levels
            for level = ADMGrid.MaxLevel:-1:2
               Prolp = obj.LevelProlongation(ADMGrid, GlobalGrids(level), level);
               % Multiply by previous objects
               ADMProlp = Prolp * ADMProlp;
            end
            
            % Last prolongation is different coz I use fine-scale ordering
            Prolp = obj.LastProlongation(ADMGrid, GlobalGrids(1));
            
            % Multiply by previous objects
            ADMProlp = Prolp * ADMProlp;
        end
        function Prolp = LevelProlongation(obj, ADMGrid, FineGrid, level)
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
             Prolp = sparse(sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nx), sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nc));
             
             % 2. Fill in top left
             Prolp(1:sum(obj.ADMmap.Nf), 1:sum(obj.ADMmap.Nf)) = speye(sum(obj.ADMmap.Nf));
             
             % 3. Fill in Bottom left
             Prolp(sum(obj.ADMmap.Nf) + 1 : end,  obj.ADMmap.Verteces) = obj.P{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexVerteces);
             
             % 4. Fill in Bottom right
             Prolp(sum(obj.ADMmap.Nf) + 1 :end, sum(obj.ADMmap.Nf) + 1 : end) = obj.P{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexNc);
        end
        function Prolp = LastProlongation(obj, ADMGrid, FineGrid)
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
             Prolp = sparse(FineGrid.N, sum(obj.ADMmap.Nf) + sum(obj.ADMmap.Nc));
             
             % 2. Fill in FS verteces of level 1
             Prolp(:,  obj.ADMmap.Verteces) = obj.P{1}(:, obj.ADMmap.OriginalIndexVerteces);
             
             % 3. Fill in coarse-scale nodes
             Prolp(:, sum(obj.ADMmap.Nf) + 1 : end) = obj.P{1}(:, obj.ADMmap.OriginalIndexNc);
             
             % 4. Fill in fine-scale nodes first
             rows = obj.ADMmap.OriginalIndexNf';
             columns = 1:sum(obj.ADMmap.Nf);
             Prolp(rows, :) = 0; % if it s fine-scale already I get rid of useless fillings
             Prolp(sub2ind(size(Prolp), rows, columns)) = 1;
        end
        function UpdateProlongationOperator(obj, FineGrid, CoarseGrid, ProductionSystem)
            % for now I do not update pressure basis functions
        end
    end
end