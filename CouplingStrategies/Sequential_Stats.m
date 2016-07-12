% Sequential coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 8 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Sequential_Stats < Coupling_Stats
    properties
        Pressure_stats
        Saturation_stats
        Pressure_timer
        Balance_timer
        Transport_timer
    end
    methods
        function obj = Sequential_Stats(MaxNTimeSteps)
           
        end
        function WriteStats(obj, itCount, n_chops)

        end
        function WriteTimers()
            
        end
    end
end