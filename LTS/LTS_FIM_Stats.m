% FIM Stats and timers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_FIM_Stats < FIM_Stats
    properties
        Complexity
    end
    methods
        function obj = LTS_FIM_Stats(MaxNTimeSteps, name)
           obj@FIM_Stats(MaxNTimeSteps, name); 
        end
        function Matrix = TimersMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.TotalTimer(1:Ndt), obj.ConstructTimer(1:Ndt), obj. SolveTimer(1:Ndt), obj.FlashTimer(1:Ndt)]';
        end
        function SaveLevelsComplexities(obj, Ndt, LTS_complexity)
            obj.Complexity(Ndt, 1:size(LTS_complexity, 2)) = LTS_complexity;
        end
    end
end