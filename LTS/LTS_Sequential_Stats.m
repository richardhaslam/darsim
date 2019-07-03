% LTS Sequential Stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_Sequential_Stats < Sequential_Stats
    properties
        NLIterLTS
        CFLGlob
        CFLLoc
        Complexity
    end
    methods
        function obj = LTS_Sequential_Stats(MaxNTimeSteps, name)
           obj@Sequential_Stats(MaxNTimeSteps, name);
           obj.NLIterLTS = zeros(MaxNTimeSteps, 1);
           obj.CFLGlob = zeros(MaxNTimeSteps, 1);
           obj.CFLLoc = zeros(MaxNTimeSteps, 1);  
           obj.NStats = 5;
           obj.Complexity = zeros(MaxNTimeSteps, 10); % I guess this should be initialized differently.
           % I just chose a random number 10 wihch is hopefully large enough
        end
        function SaveTimers(obj, Ndt, t_pressure, t_balance, t_transport)
            obj.PressureTimer(Ndt) = sum(t_pressure);
            obj.BalanceTimer(Ndt) = sum(t_balance);
            obj.TransportTimer(Ndt) = sum(t_transport);
        end
        function SaveStats(obj, Ndt, outit, NLit, CFLglob, NLitLTS, CFLloc)
             obj.OuterIter(Ndt) = outit;
             obj.NLIter(Ndt) = NLit; % this is actually the total complexity
             obj.CFLGlob(Ndt) = CFLglob;
             obj.NLIterLTS(Ndt) = NLitLTS; 
             obj.CFLLoc(Ndt) = CFLloc; 
        end
        function SaveLevelsComplexities(obj, Ndt, LTS_complexity)
            obj.Complexity(Ndt, 1:size(LTS_complexity, 2)) = LTS_complexity;
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