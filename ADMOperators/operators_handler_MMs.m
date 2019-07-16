% Operators handler for multilevel multiscale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
            obj.P = obj.ProlongationBuilders.P{1};
            for i=2:maxLevel(1)
                % Commented by Mousa obj.R = obj.R * obj.ProlongationBuilders.R{i};
                % Commented by Mousa obj.P = obj.ProlongationBuilders.P{maxLevel(1) - i + 1} * obj.P;
                obj.R = obj.ProlongationBuilders.R{i} * obj.R;
                obj.P = obj.P * obj.ProlongationBuilders.P{i};
            end
        end
    end
end