% TimeLoop driver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 7 February 20176
%Last modified: 7 February 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef end_of_sim_evaluator_totaltime < end_of_sim_evaluator
    properties
        TotalTime
        MaxNumberOfTimeSteps
    end
    methods
        function obj = end_of_sim_evaluator_totaltime(total_time, max_n_tsteps)
            obj.TotalTime = total_time;
            obj.MaxNumberOfTimeSteps = max_n_tsteps;
        end
        function End = HasSimulationEnded(obj, ProductionSystem, Time, Ndt)
            End = 0;
            if Time < obj.TotalTime || Ndt == obj.MaxNumberOfTimeSteps
            end
        end
    end
end