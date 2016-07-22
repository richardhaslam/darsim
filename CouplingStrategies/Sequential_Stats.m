% Sequential Stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 8 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Sequential_Stats < Coupling_Stats
    properties
        PressureStats
        SaturationStats
        PressureTimer
        BalanceTimer
        TransportTimer
    end
    methods
        function obj = Sequential_Stats(MaxNTimeSteps)
           obj@Coupling_Stats(MaxNTimeSteps);
           obj.PressureStats = zeros(MaxNTimeSteps, 2);
           obj.SaturationStats = zeros(MaxNTimeSteps, 2);
           obj.PressureTimer = zeros(MaxNTimeSteps, 1);
           obj.BalanceTimer = zeros(MaxNTimeSteps, 1);
           obj.TransportTimer = zeros(MaxNTimeSteps, 1);
        end
        function SaveTimers(obj)
        end
        function SaveStats(obj)
        end
        function WriteStats(obj, itCount, n_chops)
            disp('I have to write the Stats');
        end
        function WriteTimers(obj)
            disp('I have to write the timers');
            
        end
    end
end