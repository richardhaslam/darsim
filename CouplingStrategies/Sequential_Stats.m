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
        OuterIter
        NLIter
        PressureTimer
        BalanceTimer
        TransportTimer
    end
    methods
        function obj = Sequential_Stats(MaxNTimeSteps)
           obj@Coupling_Stats(MaxNTimeSteps);
           obj.OuterIter = zeros(MaxNTimeSteps, 1);
           obj.NLIter = zeros(MaxNTimeSteps, 1);
           obj.PressureTimer = zeros(MaxNTimeSteps, 1);
           obj.BalanceTimer = zeros(MaxNTimeSteps, 1);
           obj.TransportTimer = zeros(MaxNTimeSteps, 1);
           obj.NTimers = 4;
           obj.NStats = 2;
        end
        function SaveTimers(obj, Ndt, t_pressure, t_balance, t_transport)
            obj.PressureTimer(Ndt) = sum(t_pressure);
            obj.BalanceTimer(Ndt) = sum(t_balance);
            obj.TransportTimer(Ndt) = sum(t_transport);
        end
        function SaveStats(obj, Ndt, outit, NLit)
             obj.OuterIter(Ndt) = outit;
             obj.NLIter(Ndt) = NLit;
        end
         function Matrix = TimersMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.PressureTimer(1:Ndt), obj.BalanceTimer(1:Ndt), obj.TransportTimer(1:Ndt)]';
        end
        function Matrix = StatsMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.OuterIter(1:Ndt), obj.NLIter(1:Ndt)]';
        end
    end
end