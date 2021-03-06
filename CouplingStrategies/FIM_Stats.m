% FIM Stats and timers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FIM_Stats < Coupling_Stats
    properties
        NLIter
        Chops
        ConstructTimer
        SolveTimer
        FlashTimer
        CFL
    end
    methods
        function obj = FIM_Stats(MaxNTimeSteps, name)
           obj@Coupling_Stats(MaxNTimeSteps, name); 
           obj.NLIter = zeros(MaxNTimeSteps, 1);
           obj.Chops = zeros(MaxNTimeSteps, 1);
           obj.CFL = zeros(MaxNTimeSteps, 1);
           obj.ConstructTimer = zeros(MaxNTimeSteps, 1);
           obj.SolveTimer = zeros(MaxNTimeSteps, 1);
           obj.FlashTimer = zeros(MaxNTimeSteps, 1);
           obj.NTimers = 4;
           obj.NStats = 3;
        end
        function SaveStats(obj, Ndt, itCount, n_chops, cfl)
            obj.NLIter(Ndt) = itCount;
            obj.Chops(Ndt) = n_chops;
            obj.CFL(Ndt) = cfl;
        end
        function SaveTimers(obj, Ndt, t_Construct, t_Solve, t_flash)
            obj.ConstructTimer(Ndt) = sum(t_Construct);
            obj.SolveTimer(Ndt) = sum(t_Solve);
            obj.FlashTimer(Ndt) = sum(t_flash);
        end
        function Matrix = TimersMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.TotalTimer(1:Ndt), obj.ConstructTimer(1:Ndt), obj. SolveTimer(1:Ndt), obj.FlashTimer(1:Ndt)]';
        end
        function Matrix = StatsMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.Chops(1:Ndt), obj.NLIter(1:Ndt), obj.CFL(1:Ndt)]';
        end
    end
end