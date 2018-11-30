% Single Phase Coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef SinglePhase_Strategy < Coupling_Strategy
    properties
        PressureSolver
        PressureTimer
        Incompressible
        chops
    end
    methods
        function obj = SinglePhase_Strategy(name)
            obj@Coupling_Strategy(name);
        end
        function AddPressureSolver(obj, pressuresolver)
            obj.PressureSolver = pressuresolver;
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
            End = 0;
            obj.chops = 0;
            obj.Converged = 0;
            dt = obj.TimeStepSelector.ChooseTimeStep();
            obj.PressureSolver.SystemBuilder.SaveInitialState(ProductionSystem, Formulation);
            
            % Save state of current time-step (it's useful for ADM to update based on time change)
            ProductionSystem.SavePreviousState();
            
            obj.PressureSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
            obj.PressureSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
            while obj.Converged == 0  && obj.chops < 10
                if (obj.chops > 0)
                    disp('Max num. of iterations reached or stagnation detected: Time-step was chopped');
                    disp(char(5));
                    disp(['Restart Newton loop dt = ', num2str(dt)]);
                end
                % Solve pressure equation
                disp('Pressure Solver')
                disp('...............................................');
                tstart1 = tic;
                obj.PressureSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                obj.PressureTimer = toc(tstart1);
                disp('...............................................');
                if obj.PressureSolver.Converged == 0
                    dt = dt/10;
                    obj.chops = obj.chops + 1;
                    Formulation.Reset();
                    obj.PressureSolver.SystemBuilder.SetInitalGuess(ProductionSystem);
                    ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                else
                    obj.Converged = 1;
                    obj.TimeStepSelector.Update(dt, obj.PressureSolver.itCount - 1, obj.chops);
                end
                if obj.Incompressible
                    disp('Single Phase incompressible problem: No time dependency');
                    End = 1;
                    dt = 0;
                end
            end
        end
        function Summary = UpdateSummary(obj, Summary, Wells, Ndt, dt)
            Summary.CouplingStats.SaveStats(Ndt, obj.PressureSolver.itCount - 1);
            Summary.CouplingStats.SaveTimers(Ndt, obj.PressureTimer);
            Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
        end
    end
end