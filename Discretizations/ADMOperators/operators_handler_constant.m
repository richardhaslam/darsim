%  ADM operators handler base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 August 2016
%Last modified: 16 August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler_constant < operators_handler
    properties
    end
    methods
        function obj = operators_handler_constant(n)
            obj@operators_handler(n)
        end
        function BuildStaticOperators(obj, CoarseGrid, FineGrid, maxLevel, K)
            % Build Restriction and Prolongation operators for static grids
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(1), 1);
            obj.Pp{1} = obj.R{1}';
            for x = 2:maxLevel
                obj.R{x} = obj.MsRestriction(CoarseGrid(x-1), CoarseGrid(x), x);
                obj.Pp{x} = obj.R{x}';
            end
        end
        function BuildADMOperators(obj)
        end
    end
end