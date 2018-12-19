% LTS Sequential Stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_Sequential_Stats < Coupling_Stats
    properties
        OuterIter
        NLIter
        PressureTimer
        BalanceTimer
        TransportTimer
        NLIterLTS
        CFLGlob
        CFLLoc
    end
    methods
        function obj = LTS_Sequential_Stats(MaxNTimeSteps)
           obj@Coupling_Stats(MaxNTimeSteps);
           obj.OuterIter = zeros(MaxNTimeSteps, 1);
           obj.NLIter = zeros(MaxNTimeSteps, 1);
           obj.NLIterLTS = zeros(MaxNTimeSteps, 1);
           obj.CFLGlob = zeros(MaxNTimeSteps, 1);
           obj.CFLLoc = zeros(MaxNTimeSteps, 1);
           obj.PressureTimer = zeros(MaxNTimeSteps, 1);
           obj.BalanceTimer = zeros(MaxNTimeSteps, 1);
           obj.TransportTimer = zeros(MaxNTimeSteps, 1);
           obj.NTimers = 3;
           obj.NStats = 5;
        end
        function SaveTimers(obj, Ndt, t_pressure, t_balance, t_transport)
            obj.PressureTimer(Ndt) = sum(t_pressure);
            obj.BalanceTimer(Ndt) = sum(t_balance);
            obj.TransportTimer(Ndt) = sum(t_transport);
        end
        function SaveStats(obj, Ndt, outit, NLit, CFLglob, NLitLTS, CFLloc)
             obj.OuterIter(Ndt) = outit;
             obj.NLIter(Ndt) = NLit;
             obj.CFLGlob(Ndt) = CFLglob;
             obj.NLIterLTS(Ndt) = NLitLTS; 
             obj.CFLLoc(Ndt) = CFLloc; 
        end
         function Matrix = TimersMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.PressureTimer(1:Ndt), obj.BalanceTimer(1:Ndt), obj.TransportTimer(1:Ndt)]';
        end
        function Matrix = StatsMatrix(obj, Ndt)
            timesteps = 1:Ndt;
            Matrix = [timesteps', obj.OuterIter(1:Ndt), obj.NLIter(1:Ndt), obj.CFLGlob(1:Ndt), obj.NLIterLTS(1:Ndt), obj.CFLLoc(1:Ndt)]';
        end
    end
end