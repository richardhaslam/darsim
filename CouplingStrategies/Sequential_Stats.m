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
           
        end
        function WriteStats(obj, itCount, n_chops)

        end
        function WriteTimers()
            
        end
    end
end