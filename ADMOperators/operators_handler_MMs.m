% Operators handler for multilevel multiscale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 22 September 2017
%Last modified: 22 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler_MMs < operators_handler
    properties
        R
        P
        C % Correction functions
    end
    methods
        function obj = operators_handler_MMs(cf)
            obj@operators_handler(cf);
        end
        function BuildMMsOperators(obj, ProductionSystem, FluidModel, FineGrid, CrossConnections, maxLevel, CoarseGrid)
            %% 1. Compute operators for each level
            obj.ProlongationBuilders.BuildStaticOperators(ProductionSystem, FluidModel, FineGrid, CrossConnections, maxLevel, CoarseGrid);
            
            %% 2. Compute global operators (RRRRR) (PPPPP)
            obj.C = obj.ProlongationBuilders.C{1};  % Correction functions: only on level 1
            obj.R = obj.ProlongationBuilders.R{1}; 
            obj.P = obj.ProlongationBuilders.P{max(maxLevel)};
            for i=2:max(maxLevel)
                obj.R = obj.R * obj.ProlongationBuilders.R{i}; 
                obj.P = obj.ProlongationBuilders.P{max(maxLevel) - i + 1} * obj.P;
            end
        end
    end
end