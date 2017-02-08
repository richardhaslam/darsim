% TimeLoop driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 7 February 20176
%Last modified: 7 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef end_of_sim_evaluator_gascut < end_of_sim_evaluator
    properties
        Limit
    end
    methods
        function obj = end_of_sim_evaluator_gascut(limit)
            obj.Limit = limit;
        end
        function End = HasSimulationEnded(obj, ProductionSystem, Time, Ndt)
            End = 0;
            gascut = sum(ProductionSystem.Wells.Prod(1).QComponents(:,1))/...
                (sum(ProductionSystem.Wells.Prod(1).QComponents(:,1)) + sum(ProductionSystem.Wells.Prod(1).QComponents(:,2)));
            if gascut >= obj.Limit
                End = 1;
            end
        end
    end
end