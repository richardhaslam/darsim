% FIM coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 18 July 2016
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
        obj.chops = 0;
        obj.MaxChops = 10;
    end
    function dt = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
        dt = obj.TimeStepSelector.ChooseTimeStep();
        obj.NLSolver.Converged = 0;
        while (obj.NLSolver.Converged == 0 && obj.chops < obj.MaxChops) 
            % Print some info to the screen
            if (obj.chops > 0)
                disp('Maximum number of iterations was reached: time-step was chopped');
                disp(['Restart Newton loop dt = ', num2str(dt)]);
            end
            
            % Linear Solver Setup
            obj.NLSolver.LinearSolver.SetUp(ProductionSystem, DiscretizationModel);
            
            % NL solver call
            obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            
            % Chop time-step if it failed to converge
            if obj.NLSolver.Converged == 0
                dt = dt/2;
                obj.chops = obj.chops + 1;
            end
        end
        obj.TimeStepSelector.Update(dt, obj.NLSolver.itCount)
        obj.chops = 0;
    end
    function Summary = UpdateSummary(obj, Summary, Wells, Ndt, dt)
        %% Stats, timers and Injection/Production data
        Summary.CouplingStats.SaveStats(Ndt, obj.NLSolver.itCount-1, obj.chops);
        Summary.CouplingStats.SaveTimers(Ndt, obj.NLSolver.TimerConstruct, obj.NLSolver.TimerSolve, obj.NLSolver.TimerInner);
        Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
    end
end
end