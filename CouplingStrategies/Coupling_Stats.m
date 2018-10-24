% Coupling Stats base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef (Abstract) Coupling_Stats < handle
    properties 
        TotalTimer
        NTimers
        NStats
    end
    methods 
        function obj = Coupling_Stats(MaxNTimeSteps)
            obj.TotalTimer = zeros(MaxNTimeSteps, 1);
        end
        function SaveTimeStepTimer(obj, Ndt, timer)
            obj.TotalTimer(Ndt) = timer; 
        end
    end
    methods (Abstract)
        obj = SaveStats(obj);
        obj = SaveTimers(obj);
    end
end