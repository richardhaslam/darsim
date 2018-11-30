% Sequential Stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef SinglePhase_Stats < Coupling_Stats
    properties
        NLIter
        PressureTimer
    end
    methods
        function obj = SinglePhase_Stats(MaxNTimeSteps)
           obj@Coupling_Stats(MaxNTimeSteps);
           obj.NLIter = zeros(MaxNTimeSteps, 1);
           obj.PressureTimer = zeros(MaxNTimeSteps, 1);
           obj.NTimers = 1;
           obj.NStats = 1;
        end
        function SaveTimers(obj, Ndt, t_pressure)
            obj.PressureTimer(Ndt) = sum(t_pressure);
        end
        function SaveStats(obj, Ndt, NLit)
             obj.NLIter(Ndt) = NLit;
        end
         function Matrix = TimersMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.PressureTimer(1:Ndt)]';
        end
        function Matrix = StatsMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.NLIter(1:Ndt)]';
        end
    end
end