%  ADM operators handler constant
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
        function obj = operators_handler_constant(n, cf)
            obj@operators_handler(n, cf)
        end
        function BuildStaticOperators(obj, CoarseGrid, FineGrid, maxLevel, K, S, FluidModel)
            % Build Restriction and Prolongation operators for static grids
            disp('Building Restriction 1 ');
            obj.R{1} = obj.MsRestriction(FineGrid, CoarseGrid(1));
            disp('Building Prolongation 1 ');
            obj.Pp{1} = obj.R{1}';
            for x = 2:maxLevel
                disp(['Building Restriction ', num2str(x)]);
                obj.R{x} = obj.MsRestriction(CoarseGrid(x-1), CoarseGrid(x));
                disp(['Building Pronlongation ', num2str(x)]);
                obj.Pp{x} = obj.R{x}';
            end
        end
        function ADMProlongation(obj, ADMGrid, FineGrid, CoarseGridid)
            % Since it s constant interpolation it is just transpose(R)
            obj.ADMProlp = obj.ADMRest';
            obj.ADMProls = obj.ADMRest';
        end
    end
end