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
        % Phase Mobilities and total Mobility
        Formulation.ComputeTotalMobility(ProductionSystem, FluidModel);
        
        % Save state of current time-step
        ProductionSystem.SavePreviousState();
        % Set up linear solver
        while obj.Converged == 0 && obj.Chops < 10
            obj.itCount = 1;
            obj.PressureSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
            if (obj.Chops == 0)
                obj.PressureSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
                obj.TransportSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
            end
            % Outer loop (Flow transport coupling)
            while obj.Converged == 0 && obj.itCount < obj.MaxIter
                % copy state to check outer convergence
                State_old = status();
                State_old.CopyProperties(ProductionSystem.Reservoir.State);
                disp(newline);
                disp(['Outer iteration: ', num2str(obj.itCount)]);
                
                %% 1. Solve pressure equation
                disp('...............................................');
                disp('Pressure Solver')
                disp('...............................................');
                tstart1 = tic;
                % Save initial State
                obj.PressureSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
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
                    dt_sub = obj.TimeStepSelector.StableTimeStep(ProductionSystem, DiscretizationModel, FluidModel, Formulation.Utot);
                end
                n_subT = max(floor(dt/dt_sub), 1);
                dt_sub = dt/n_subT;
                
                disp('Transport Solver');
                disp('...............................................');
                tstart3 = tic;
                % 3.2 Solve Transport until dt of pressure is reached
                Sync = 0;
                SumDt_sub = 0;
                SubT_counter = 1;
                Temp = status();
                Temp.CopyProperties(ProductionSystem.Reservoir.State_old);
                chops = 0;
                while Sync == 0
                    % Save new old state
                    disp(['SubTime-step ', num2str(SubT_counter)]);
                    obj.TransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt_sub);
                    obj.TransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt_sub);
                    obj.NLiter = obj.NLiter + obj.TransportSolver.itCount - 1;
                    if obj.TransportSolver.Converged == 0
                        disp(['Cutting subtimestep ', num2str(SubT_counter)]);
                        % Reset initial guess to what it was before
                        obj.TransportSolver.SystemBuilder.SetInitalGuess(ProductionSystem);
                        dt_sub = dt_sub / 10;
                        chops = chops + 1;
                        if chops >= 5
                            break
                        end
                    else
                        ProductionSystem.SavePreviousState();
                        SubT_counter = SubT_counter + 1;
                        SumDt_sub = SumDt_sub + dt_sub; 
                        if norm(SumDt_sub - dt)/dt < 1e-3 
                            Sync = 1;
                            ProductionSystem.Reservoir.State_old.CopyProperties(Temp);
                        end
                    end
                end
                obj.TransportTimer(obj.itCount) = toc(tstart3);
                disp('...............................................');
                
                %% 4. Check outer-loop convergence
                if Sync == 1
                    obj.Converged = obj.ConvergenceChecker.Check(ProductionSystem.Reservoir.State, State_old);
                    obj.itCount = obj.itCount + 1;
                else
                    break
                end
            end
            
            if obj.Converged == 0
                obj.Chops = obj.Chops + 1; 
                % Reset Initial guess
                ProductionSystem.Reservoir.State.CopyProperties(ProductionSystem.Reservoir.State_old);
                % UpdateWells
                ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                disp(['Chopping time-step and restarting outer-loop with dt = ', num2str(dt)]);
                dt = dt/2;
                obj.TransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
            end
        end
        % Compute phase fluxes after converged solution
        ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        obj.TimeStepSelector.Update(dt, obj.itCount - 1, obj.Chops);
    end
    function UpdateSummary(obj, Summary, Wells, Ndt, dt)
        Summary.CouplingStats.SaveStats(Ndt, obj.itCount - 1, obj.NLiter);
        Summary.CouplingStats.SaveTimers(Ndt, obj.PressureTimer, obj.BalanceTimer, obj.TransportTimer);
        Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
    end
end
end