%  ADM operators handler MS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 6 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler_MS < operators_handler
    properties
        BFUpdater
    end
    methods
        function obj = operators_handler_MS(n, cf)
            obj@operators_handler(n, cf)
        end
        function BuildStaticOperators(obj, CoarseGrid, FineGrid, maxLevel, K, s, FluidModel)
            % Remove high contrast to avoid spikes
            lambdaMax = max(K(:,1));
            K(K(:,1)./lambdaMax < 10^-1, 1) = 10^-1*lambdaMax;
            K(K(:,2)./lambdaMax < 10^-1, 2) = 10^-1*lambdaMax;
            % Build Pressure system
            obj.BFUpdater.ConstructPressureSystem(FineGrid, K, s, FluidModel);
            %Build static restriction operator (FV)
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(1));
            % Build Prolongation operator
            obj.Pp{1} = obj.BFUpdater.MsProlongation(FineGrid, CoarseGrid(1), CoarseGrid(1).CoarseFactor);
            %Build first coarse system (with MsFE)
            obj.BFUpdater.A = obj.Pp{1}' * obj.BFUpdater.A * obj.Pp{1};
            obj.BFUpdater.TransformIntoTPFA(CoarseGrid(1).Nx);
            for x = 2:maxLevel
                % Build static restriction operator (FV)
                obj.R{x} = obj.MsRestriction(CoarseGrid(x-1), CoarseGrid(x));
                % Build Prolongation operator
                obj.Pp{x} = obj.BFUpdater.MsProlongation(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x).CoarseFactor./CoarseGrid(x-1).CoarseFactor);
                %Build coarse system (with MsFE)
                obj.BFUpdater.A = obj.Pp{x}' * obj.BFUpdater.A * obj.Pp{x};
                obj.BFUpdater.TransformIntoTPFA(CoarseGrid(x).Nx);
            end
        end
        function ADMProlongation(obj, ADMGrid, FineGrid, CoarseGrid)
            % Pressure prolongation
            obj.ADMProlp = 1;
            
            % Loop over the levels
            for level = ADMGrid.MaxLevel:-1:2
               Prolp = obj.LevelProlongation(ADMGrid, CoarseGrid(level-1), level);
               % Multiply by previous objects
               obj.ADMProlp = Prolp * obj.ADMProlp;
            end
            
            % Last prolongation is different coz I use fine-scale ordering
            Prolp = obj.LastProlongation(ADMGrid, FineGrid);
            
            % Multiply by previous objects
            obj.ADMProlp = Prolp * obj.ADMProlp;
            
            % Saturation prolongation: transpose(R)
            obj.ADMProls = obj.ADMRest'; 
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
             Prolp = zeros(obj.ADMmap.Nf + obj.ADMmap.Nx, obj.ADMmap.Nf + obj.ADMmap.Nc);
             
             % 2. Fill in top left
             Prolp(1:obj.ADMmap.Nf, 1:obj.ADMmap.Nf) = eye(obj.ADMmap.Nf);
             
             % 3. Fill in Bottom left
             Prolp(obj.ADMmap.Nf + 1 : end,  obj.ADMmap.Verteces) = obj.Pp{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexVerteces);
             
             % 4. Fill in Bottom right
             Prolp(obj.ADMmap.Nf + 1 :end, obj.ADMmap.Nf + 1 : end) = obj.Pp{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexNc);
             
             % 5. Make it sparse
             Prolp = sparse(Prolp);
        end
        function Prolp = LastProlongation(obj, ADMGrid, FineGrid, CoarseGrid)
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % For a the firse level level the prolongation operator looks like this             
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
             Prolp = zeros(FineGrid.N, obj.ADMmap.Nf + obj.ADMmap.Nc);
             
             % 2. Fill in FS verteces of level 1
             Prolp(:,  obj.ADMmap.Verteces) = obj.Pp{1}(:, obj.ADMmap.OriginalIndexVerteces);
             
             % 3. Fill in coarse-scale nodes
             Prolp(:, obj.ADMmap.Nf + 1 : end) = obj.Pp{1}(:, obj.ADMmap.OriginalIndexNc);
             
             % 4. Fill in fine-scale nodes first
             rows = obj.ADMmap.OriginalIndexNf';
             columns = 1:obj.ADMmap.Nf;
             Prolp(rows,:) = 0; % if it s fine-scale already I get rid of useless fillings
             Prolp(sub2ind(size(Prolp), rows, columns)) = 1;
             
             % 5. Make it sparse
             Prolp = sparse(Prolp);
        end
    end
end