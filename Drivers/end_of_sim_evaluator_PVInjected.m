% TimeLoop driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 6 April 2017
%Last modified: 6 April 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef end_of_sim_evaluator_PVInjected < end_of_sim_evaluator
    properties
        Limit
    end
    methods
        function obj = end_of_sim_evaluator_PVInjected(total_time, max_n_tsteps, limit)
            obj@end_of_sim_evaluator(total_time, max_n_tsteps);
            obj.Limit = limit;
        end
        function End = HasSimulationEnded(obj, End, Summary, ProductionSystem, Time, Ndt)
            % It has to be done better. For now it whould be fine like
            % this. In theory it would be better to store volumetric
            % injection in the wells data.
            PVinjected = sum(Summary.WellsData.Injection.Phases(Ndt, :, 1)) / (max(ProductionSystem.Wells.Inj(1).rho(:,1)) * ProductionSystem.Reservoir.TotalPV);
            if PVinjected >= obj.Limit || Time >= obj.TotalTime || Ndt == obj.MaxNumberOfTimeSteps
                End = 1;
            end
        end
    end
end