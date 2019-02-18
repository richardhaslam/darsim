% Sequential Stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Sequential_Stats < Coupling_Stats
    properties
        OuterIter
        NLIter
        CFLVal
        PressureTimer
        BalanceTimer
        TransportTimer
    end
    methods
        function obj = Sequential_Stats(MaxNTimeSteps, name)
           obj@Coupling_Stats(MaxNTimeSteps, name);
           obj.OuterIter = zeros(MaxNTimeSteps, 1);
           obj.NLIter = zeros(MaxNTimeSteps, 1);
           obj.CFLVal = zeros(MaxNTimeSteps, 1);
           obj.PressureTimer = zeros(MaxNTimeSteps, 1);
           obj.BalanceTimer = zeros(MaxNTimeSteps, 1);
           obj.TransportTimer = zeros(MaxNTimeSteps, 1);
           obj.NTimers = 2;
           obj.NStats = 2;
        end
        function SaveTimers(obj, Ndt, t_pressure, t_balance, t_transport)
            obj.PressureTimer(Ndt) = sum(t_pressure);
            obj.BalanceTimer(Ndt) = sum(t_balance);
            obj.TransportTimer(Ndt) = sum(t_transport);
        end
        function SaveStats(obj, Ndt, outit, NLit, CFLval)
             obj.OuterIter(Ndt) = outit;
             obj.NLIter(Ndt) = NLit;
             obj.CFLVal(Ndt) = CFLval;
        end
         function Matrix = TimersMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.PressureTimer(1:Ndt), obj.TransportTimer(1:Ndt)]';
        end
        function Matrix = StatsMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.OuterIter(1:Ndt), obj.NLIter(1:Ndt)]';
        end
    end
end