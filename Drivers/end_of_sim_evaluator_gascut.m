% TimeLoop driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef end_of_sim_evaluator_gascut < end_of_sim_evaluator
    properties
        Limit
    end
    methods
        function obj = end_of_sim_evaluator_gascut(total_time, max_n_tsteps, limit)
            obj@end_of_sim_evaluator(total_time, max_n_tsteps);
            obj.Limit = limit;
            obj.Limit = limit;
        end
        function End = HasSimulationEnded(obj, End, Summary, ProductionSystem, Time, Ndt)
            gascut = sum(Summary.WellsData.Production.Phases(Ndt, :, 1))/sum(Summary.WellsData.Injection.Phases(Ndt, :, 1));
            if gascut >= obj.Limit
                End = 1;
            end
        end
    end
end