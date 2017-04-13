% Single Phase Coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 2 March 2017
%Last modified: 2 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef SinglePhase_Strategy < Coupling_Strategy
properties
    PressureSolver
    PressureTimer
end
methods
    function obj = SinglePhase_Strategy(name)
        obj@Coupling_Strategy(name);
    end
    function AddPressureSolver(obj, pressuresolver)
        obj.PressureSolver = pressuresolver;
    end
    function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)        
        dt = obj.TimeStepSelector.ChooseTimeStep();
        obj.PressureSolver.SystemBuilder.SaveInitialState(ProductionSystem, Formulation);
        % Solve pressure equation
        disp('Pressure Solver')
        disp('...............................................');
        tstart1 = tic;
        obj.PressureSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
        obj.PressureTimer = toc(tstart1);
        disp('...............................................');
        dt = 0;
        End = 1;
        obj.Converged = 1;
    end
    function Summary = UpdateSummary(obj, Summary, Wells, Ndt, dt)
        Summary.CouplingStats.SaveStats(Ndt, obj.PressureSolver.itCount - 1);
        Summary.CouplingStats.SaveTimers(Ndt, obj.PressureTimer);
        Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
    end
end
end