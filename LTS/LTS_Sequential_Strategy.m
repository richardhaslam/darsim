% LTS SEQUENTIAL strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author:Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_Sequential_Strategy < Sequential_Strategy
    properties
        LTSTransportSolver
        RefCellsSelector
        LTSTransportTimer
        NofRef = 4;
    end
    methods
        function AddLTSTransportSolver(obj, LTStransportsolver)
            obj.LTSTransportSolver = LTStransportsolver;
        end
        function obj = LTS_Sequential_Strategy(name)
            obj@Sequential_Strategy(name);
            obj.RefCellsSelector = RefCellsSelector();
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
                    obj.LTSTransportSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);

                end
                % Outer loop (Flow transport coupling)
                while obj.Converged == 0 && obj.itCount < obj.MaxIter
                    % copy state to check outer convergence
                    State_old = status();
                    State_old.CopyProperties(ProductionSystem.Reservoir.State);
                    disp(newline);
                    Formulation.ComputeTotalFluxes(ProductionSystem, DiscretizationModel);
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
                    
                    disp('Transport Solver');
                    disp('...............................................');
                    tstart3 = tic;
                    % 3.2 Solve Transport tentative global step
       
                    obj.TransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
                    obj.TransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                    if obj.TransportSolver.Converged == 0
                        error('DARSim2 error: transport not converged')
                    end
                    obj.TransportTimer(obj.itCount) = toc(tstart3);
                    disp('...............................................');
                    
                    % 3.3 Compute active cells
                   
                    obj.RefCellsSelector.SelectRefCells(ProductionSystem, DiscretizationModel.ReservoirGrid, Formulation);
                    
                    
                    %if there is at least one active cell we need to
                    %recompute the new smaller system 
                    if sum(sum(obj.RefCellsSelector.ActCells)) > 0

                        DiscretizationModel.ReservoirGrid.ActiveTime = obj.RefCellsSelector.ActCells;
                        State_global = status(); 
                        State_global.CopyProperties(ProductionSystem.Reservoir.State);
                        dtRef = dt / obj.NofRef;
                        % Compute the numerical fluxes used as boundary
                        % values between the accepted and rejected area.
                        obj.RefCellsSelector.ComputeBoundaryValues(DiscretizationModel, Formulation);
                        %Set initial values for the saturation
                        obj.TransportSolver.SystemBuilder.SetInitalGuess(ProductionSystem);
                        
                        % we sum up all the timer for the refinemets
                        obj.LTSTransportTimer(obj.itCount) = 0;
                        disp('Transport Solver Sub-rebinements');
                        disp('...............................................');
                        for itSub = 1 : obj.NofRef
                            % Solve the transport eq for the subref cells
                            
                            disp(['SubRef step: ', num2str(itSub)]);
                            disp('...............................................'); 
                            tstart4 = tic;
                            
                            obj.LTSTransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dtRef, obj.RefCellsSelector);
                            obj.LTSTransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dtRef, obj.RefCellsSelector);
                            if obj.LTSTransportSolver.Converged == 0
                                disp('Sub ref not converged')
                                break;
                            end
                            
                            ProductionSystem.SavePreviousState();
                            %for each subref    
                            obj.LTSTransportTimer(obj.itCount) = obj.LTSTransportTimer(obj.itCount) + toc(tstart4);
                            disp('...............................................');
                            
                        end
                        obj.LTSTransportSolver.SynchronizeProperties(ProductionSystem, State_global, obj.RefCellsSelector);  
                    else
                        obj.LTSTransportSolver.Converged = 1;
                        obj.LTSTransportTimer(obj.itCount) = 0;
                        DiscretizationModel.ReservoirGrid.ActiveTime = ones(DiscretizationModel.ReservoirGrid.N,1);
                    end
                    obj.Converged = obj.ConvergenceChecker.Check(ProductionSystem.Reservoir.State, State_old);
                    obj.itCount = obj.itCount + 1;
                    if (obj.LTSTransportSolver.Converged == 0 || obj.LTSTransportSolver.Converged == 0 || obj.PressureSolver.Converged == 0)
                        obj.Converged = 0;
                        
                        obj.Chops = obj.Chops + 1;
                        % Reset Initial guess
                        ProductionSystem.Reservoir.State.CopyProperties(ProductionSystem.Reservoir.State_old);
                        % UpdateWells
                        ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                        disp(['Chopping time-step and restarting outer-loop with dt = ', num2str(dt)]);
                        dt = dt/2;
                        obj.TransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
                        break;
                    end
                end 
                %ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                ProductionSystem.SavePreviousState();
                obj.TimeStepSelector.Update(dt, obj.itCount - 1, obj.Chops);
            end
        end
    end
end