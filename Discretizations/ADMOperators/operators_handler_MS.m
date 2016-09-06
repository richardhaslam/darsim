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
        function obj = operators_handler_MS(n)
            obj@operators_handler(n)
        end
        function BuildStaticOperators(obj, CoarseGrid, FineGrid, maxLevel, K, s, FluidModel)
            obj.BFUpdater.ConstructPressureSystem(FineGrid, K, s, FluidModel);
            %Build static restriction operator (FV)
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(1));
            % Build Prolongation operator
            obj.Pp{1} = obj.BFUpdater.MsProlongation(FineGrid, CoarseGrid(1), CoarseGrid(1).CoarseFactor);
            %Build first coarse system (with MsFE)
            obj.BFUpdater.A = obj.Pp{1}' * obj.BFUpdater.A * obj.Pp{1};
            for x = 2:maxLevel
                % Build static restriction operator (FV)
                obj.R{x} = obj.MsRestriction(CoarseGrid(x-1), CoarseGrid(x));
                % Build Prolongation operator
                obj.Pp{x} = obj.BFUpdater.MsProlongation(CoarseGrid(x-1), CoarseGrid(x), CoarseGrid(x).CoarseFactor./CoarseGrid(x-1).CoarseFactor);
                %Build coarse system (with MsFE)
                obj.BFUpdater.A = obj.Pp{x}' * obj.BFUpdater.A * obj.Pp{x};
            end
        end
        function ADMProlongation(obj)
            % Pressure prolongation
            obj.ADMProlp = 1;
            
            % Loop over the levels
            for level = ADMGrid.MaxLevel:-1:1
               Prolp = LevelProlongation();
               % Multiply by previous objects
               obj.ADMProlp = Prolp * obj.ADMProlp;
            end
            
            
            
            % Saturation prolongation: transpose(R)
            obj.ADMProls = obj.ADMRest'; 
        end
        function Prolp = LevelProlongation(obj, level)
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
             obj.ADMmap.update();
             
             % 1. Build the object
             Prolp = zeros(obj.ADMmap.Nf + obj.ADMmap.Nx, obj.ADMmap.Nf + obj.ADMmap.Nc);
             
             % 2. Fill in top left
             Prolp(Nf, Nf) = eye(Nf);
             
             % 3. Fill in Bottom left
             Prolp(obj.ADMmap.Nf + 1 : end,  obj.ADMmap.Verteces) = obj.Pp{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexVerteces);
             
             % 4. Fill in Bottom right
             Prolp(obj.ADMmap.Nf + 1 :end, obj.ADMmap.Nf + 1: end) = obj.Pp{level}(obj.ADMmap.OriginalIndexNx, obj.ADMmap.OriginalIndexNc);
             
             % 5. Make it sparse
             Prolp = sparse(Prolp);
        end
    end
end