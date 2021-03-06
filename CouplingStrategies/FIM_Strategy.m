% FIM coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FIM_Strategy < Coupling_Strategy
    properties
        NLSolver
        chops
        MaxChops
        CFL
    end
    methods
        function obj = FIM_Strategy(name, NONLinearSolver)
            obj@Coupling_Strategy(name);
            obj.NLSolver = NONLinearSolver;
            obj.chops = 0;
            obj.MaxChops = 15;
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
            
            % Initialise
            dt = obj.TimeStepSelector.ChooseTimeStep();
            obj.Converged = 0;
            obj.chops = 0;
            End = 0;
            
            DT = obj.TimeStepSelector.Sec2DHMS(dt);
            disp(['dt= ' num2str(dt) ' sec (', num2str(DT.Days), ' days : ', num2str(DT.Hours), ' hrs : ', num2str(DT.Minutes), ' mins : ', num2str(DT.Seconds), ' sec)']);
            
            % Set Up non-linear solver
            obj.NLSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
            obj.NLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
            
            % Save state of current time-step (it's useful for ADM to update based on time change)
            ProductionSystem.SavePreviousState();
            while (obj.Converged == 0 && obj.chops < obj.MaxChops)
                % Print some info to the screen
                if (obj.chops > 0)
                    disp('Max num. of iterations reached or stagnation detected: Time-step was chopped');
                    disp(newline);
                    disp(['Restart Newton loop dt = ', num2str(dt)]);
                end
                
                % NL solver call
                obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                
                % Chop time-step if it failed to converge
                if obj.NLSolver.Converged == 0
                    dt = dt/2;
                    obj.chops = obj.chops + 1;
                    Formulation.Reset();
                    obj.NLSolver.SystemBuilder.SetInitalGuess(ProductionSystem);
                    ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                    obj.NLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
                end
                obj.Converged = obj.NLSolver.Converged;
            end
            Formulation.ComputeTotalFluxes(ProductionSystem, DiscretizationModel);
            obj.CFL = Formulation.ComputeCFLNumber(ProductionSystem, DiscretizationModel, dt);
            % disp(['CFL = ', num2str(obj.CFL)]);
            obj.TimeStepSelector.Update(dt, obj.NLSolver.itCount - 1, obj.chops);
        end
        function Summary = UpdateSummary(obj, Summary, Wells, Ndt, dt)
            %% Stats, timers and Injection/Production data
            Summary.CouplingStats.SaveStats(Ndt, obj.NLSolver.itCount-1, obj.chops, obj.CFL);
            Summary.CouplingStats.SaveTimers(Ndt, obj.NLSolver.TimerConstruct, obj.NLSolver.TimerSolve, obj.NLSolver.TimerInner);
            Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
        end
    end
end