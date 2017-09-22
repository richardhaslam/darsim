% Sequential coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 23 July 2016
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
        obj.itCount = 1;
        obj.Converged = 0;
        obj.NLiter = 0;
        
        % Timers
        obj.PressureTimer = zeros(obj.MaxIter,1);
        obj.BalanceTimer = zeros(obj.MaxIter, 1);
        obj.TransportTimer = zeros(obj.MaxIter, 1);
        
        % Pressure timestep
        dt = obj.TimeStepSelector.ChooseTimeStep();
        % Phase Mobilities and total Mobility
        Formulation.ComputeTotalMobility(ProductionSystem, FluidModel);
        
        % Choose Grid resolution for this time-step (does nothing for FS)
        DiscretizationModel.SelectADMGrid(ProductionSystem);
        
        % Save initial State
        obj.PressureSolver.SystemBuilder.SaveInitialState(ProductionSystem, Formulation);
        obj.TransportSolver.SystemBuilder.SaveInitialState(ProductionSystem, Formulation);
        % Save state of current time-step (it's useful for ADM to update based on time change)
        ProductionSystem.SavePreviousState();
        while obj.Converged == 0 && obj.itCount < obj.MaxIter
            % copy state to check outer convergence
            State_old = status();
            State_old.CopyProperties(ProductionSystem.Reservoir.State);
            disp(newline);
            disp(['Outer iteration: ', num2str(obj.itCount)]);
            %% 1. Solve pressure equation
            disp('Pressure Solver')
            disp('...............................................');
            tstart1 = tic;
            obj.PressureSolver.LinearSolver.SetUp(DiscretizationModel);
            obj.PressureSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            obj.PressureTimer(obj.itCount) = toc(tstart1);
            disp('...............................................');
            %% 2. Compute total fluxes 
            Formulation.UpWindAndPhaseRockFluxes(DiscretizationModel, FluidModel.Phases, ProductionSystem);
            Formulation.ComputeTotalFluxes(ProductionSystem, DiscretizationModel);
            % Check that velocity field is conservative
            tstart2 = tic;
            conservative = Formulation.CheckMassConservation(DiscretizationModel.ReservoirGrid);
            if ~conservative
                error('DARSim2 error: mass balance not respected');
            end
            obj.BalanceTimer(obj.itCount) = toc(tstart2);
             
            %% 3. Solve transport            
            % 3.1 Choose stable timestep
            if obj.itCount == 1
                dt = obj.TimeStepSelector.StableTimeStep(ProductionSystem, DiscretizationModel, FluidModel, Formulation.Utot);
            end
            disp('Transport Solver');
            disp('...............................................');
            tstart3 = tic;
            obj.TransportSolver.Converged = 0;
            % 3.2 Solve
            TempState = status();
            TempState.CopyProperties(ProductionSystem.Reservoir.State);
            while obj.TransportSolver.Converged == 0
                obj.TransportSolver.LinearSolver.SetUp(DiscretizationModel);
                obj.TransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                obj.NLiter = obj.NLiter + obj.TransportSolver.itCount - 1;
                if obj.TransportSolver.Converged == 0
                    ProductionSystem.Reservoir.State.CopyProperties(TempState);
                    % UpdateWells
                    ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                    disp(['Chopping time-step and restarting Newton loop with dt = ', num2str(dt)])
                    dt = dt/10;
                end
            end
            obj.TransportTimer(obj.itCount) = toc(tstart3);
            disp('...............................................');
            
            %% 4. Check outer-loop convergence
            obj.Converged = obj.ConvergenceChecker.Check(ProductionSystem.Reservoir.State, State_old);
            
            obj.itCount = obj.itCount + 1;
        end
        % Compute phase fluxes after converged solution
        ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
    end
    function UpdateSummary(obj, Summary, Wells, Ndt, dt)
        Summary.CouplingStats.SaveStats(Ndt, obj.itCount - 1, obj.NLiter);
        Summary.CouplingStats.SaveTimers(Ndt, obj.PressureTimer, obj.BalanceTimer, obj.TransportTimer);
        Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
    end
end
end