% Sequential coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Sequential_Strategy < Coupling_Strategy
    properties
        PressureSolver
        TransportSolver
        MaxIter
        itCount
        ConvergenceChecker
        NLiter
        PressureTimer
        BalanceTimer
        TransportTimer
        Chops
    end
    methods
        function obj = Sequential_Strategy(name)
            obj@Coupling_Strategy(name);
        end
        function AddPressureSolver(obj, pressuresolver)
            obj.PressureSolver = pressuresolver;
        end
        function AddTransportSolver(obj, transportsolver)
            obj.TransportSolver = transportsolver;
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
            End = 0;
            
            % This is the outer loop
            obj.Converged = 0;
            obj.NLiter = 0;
            obj.Chops = 0;
            
            % Timers
            obj.PressureTimer = zeros(obj.MaxIter,1);
            obj.BalanceTimer = zeros(obj.MaxIter, 1);
            obj.TransportTimer = zeros(obj.MaxIter, 1);
            
            % Pressure timestep
            dt = obj.TimeStepSelector.ChooseTimeStep();
            
            % Save state of current time-step
            ProductionSystem.SavePreviousState();
            % Set up linear solver
            while obj.Converged == 0 && obj.Chops < 10
                obj.itCount = 1;
                % Save initial State (only at first outerloop iteration)
                obj.PressureSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
                
                if (obj.Chops == 0)
                    % This is only for ADM
                    obj.PressureSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
                    obj.TransportSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
                end
                
                % Outer loop (Flow-transport coupling)
                while obj.Converged == 0 && obj.itCount <= obj.MaxIter
                    % copy state to check outer convergence
                    State_old = status();
                    State_old.CopyProperties(ProductionSystem.Reservoir.State);
                    disp(['Outer iteration: ', num2str(obj.itCount)]);
                    
                    %% 1. Solve pressure equation
                    disp('...............................................');
                    disp('Pressure Solver')
                    disp('...............................................');
                    tstart1 = tic;
                    obj.PressureSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                    obj.PressureTimer(obj.itCount) = toc(tstart1);
                    disp('...............................................');
                    
                    %% 2. Compute total fluxes
                    Formulation.UpWindAndPhaseRockFluxes(DiscretizationModel, FluidModel.Phases, ProductionSystem);
                    Formulation.ComputeTotalFluxes(ProductionSystem, DiscretizationModel);
                    % Check that velocity field is conservative
                    tstart2 = tic;
%                     % conservative = Formulation.CheckMassConservation(DiscretizationModel.ReservoirGrid);
%                     if ~conservative
%                         error('DARSim2 error: mass balance not respected');
%                     end
                    obj.BalanceTimer(obj.itCount) = toc(tstart2);
                    
                    %% 3. Solve transport
                    % 3.1 Choose stable timestep
                    if obj.itCount == 1 && obj.Chops == 0
                        dt = obj.TimeStepSelector.StableTimeStep(ProductionSystem, DiscretizationModel, FluidModel, Formulation.Utot);
                    end
                    % 3.2 Solve transport
                    disp('Transport Solver');
                    disp('...............................................');
                    tstart3 = tic;
                    if obj.itCount == 1
                        % Save new old state
                        obj.TransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
                    end
                    obj.TransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                    obj.NLiter = obj.NLiter + obj.TransportSolver.itCount - 1;
                    obj.TransportTimer(obj.itCount) = toc(tstart3);
                    disp('...............................................');
                    disp(newline);
                    if obj.TransportSolver.Converged == 0
                        disp('Transport solver failed to converge!');
                        obj.itCount = obj.MaxIter+1;
                    else
                        %% 4. Check outer-loop convergence
                        obj.Converged = obj.ConvergenceChecker.Check(ProductionSystem.Reservoir.State, State_old);
                        obj.itCount = obj.itCount + 1;
                    end
                end
                
                % If it did not converge we chop the time-step and try
                % again
                if obj.Converged == 0
                    dt = dt/2;
                    obj.Chops = obj.Chops + 1;
                    % Reset Initial guess
                    ProductionSystem.Reservoir.State.CopyProperties(ProductionSystem.Reservoir.State_old);
                    % UpdateWells
                    ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                    disp(['Chopping time-step and restarting outer-loop with dt = ', num2str(dt)]);
                    disp(newline);
                end
            end
            % Compute phase fluxes after converged solution
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            obj.TimeStepSelector.UpdateSequential(dt, obj.itCount - 1, obj.Chops);
        end
        function UpdateSummary(obj, Summary, Wells, Ndt, dt)
            Summary.CouplingStats.SaveStats(Ndt, obj.itCount - 1, obj.NLiter);
            Summary.CouplingStats.SaveTimers(Ndt, obj.PressureTimer, obj.BalanceTimer, obj.TransportTimer);
            Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
        end
    end
end