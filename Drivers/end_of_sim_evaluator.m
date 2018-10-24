% EndOfSimEvaluator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef end_of_sim_evaluator 
    properties
        TotalTime
        MaxNumberOfTimeSteps
    end
    methods
        function obj = end_of_sim_evaluator(total_time, max_n_tsteps)
            obj.TotalTime = total_time;
            obj.MaxNumberOfTimeSteps = max_n_tsteps;
        end
        function End = HasSimulationEnded(obj, End, Summary, ProductionSystem, Time, Ndt)
            Ndt = Ndt - 1;
            if Time >= obj.TotalTime || Ndt == obj.MaxNumberOfTimeSteps
                End = 1;
            end
        end
    end
end

