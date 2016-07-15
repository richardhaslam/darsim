% FIM coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 8 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FIM_Strategy < Coupling_Strategy
properties
    NLSolver
    TimeStepSelector
    chops
    MaxChops
end
methods
    function obj = FIM_Strategy(name, NONLinearSolver)
        obj@Coupling_Strategy(name);
        obj.NLSolver = NONLinearSolver;
    end
    function [ProductionSystem, dt] = SolveTimeStep(obj, ProductionSystem, DiscretizationModel, Formulation, maxDt)
        dt = obj.TimeStepSelector(maxDt);
        while (obj.NLSolver.Converged == 0 && obj.chops < obj.MaxChops) 
            % Print some info to the screen
            if (obj.chops > 0)
                disp('Maximum number of iterations was reached: time-step was chopped');
                disp(['Restart Newton loop dt = ', num2str(dt)]);
            end
            
            % Linear Solver Setup
            obj.NLSolver.LinearSolver.setup(ProductionSystem, DiscretizationModel);
            
            % NL solver call
            [ProductionSystem] = obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            
            % Chop time-step if it failed to converge
            if obj.NLSolver.Converged == 0
                dt = dt/2;
                obj.chops = obj.chops + 1;
            end
        end
        obj.TimeStepSelector.Update(dt, obj.NLSolver.itCount)
    end
    function Summary = UpdateSummaryStats(Summary, Ndt)
        %% Stats, timers and Injection/Production data
        Summary.CouplingStats.SaveStats(obj.NLSolver.Ndt, obj.NLSolver.itCount-1, obj.chops);
        Summary.CouplingStats.SaveTimers(Ndt, obj.NLSolver.TimerConstruct, obj.NLSolver.TimerSolve, obj.NLSolver.TimerInner);
        Summary.SaveWellsData(Ndt+1, ProductionSystem.Wells.Inj, ProductionSystem.Wells.Prod, dt);
    end
end
end