%  ADM operators handler base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
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
        function ADMProlongation(obj, level)
            % Pressure prolongation
            Prolp = obj.ADMRest{level}';
            % Copy MsP in ADM prolongation operator
            Prolp(Nf:end, ADMmap.Verteces(:,1)) = obj.Pp{level}(ADMmap.OriginalIndex(Nf:end), ADMmap.Verteces(:,2));
            Prolp(Nf:end, Nf + Nc) = obj.Pp{level}(ADMmap.OriginalIndex(Nf:end), ADMmap.OriginalIndex(Nf + Nc:end));
            
            % Copy it in the right object
            obj.ADMProlp{level} = Prolp;
            % Saturation prolongation: transpose(R)
            obj.ADMProls{level} = obj.ADMRest{level}'; 
        end
    end
end