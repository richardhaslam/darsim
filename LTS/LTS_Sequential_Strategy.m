% LTS Sequential strategy
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
        NLiterLTS
        CFLGlobal
        CFLLocal
        ActCellsSummary
        StatesSummary
        NofRef = 4;
    end
    methods 
        function obj = LTS_Sequential_Strategy(name,tol)
            obj@Sequential_Strategy(name);
            obj.RefCellsSelector = RefCellsSelector(tol);
        end
        function AddLTSTransportSolver(obj, LTStransportsolver)
            obj.LTSTransportSolver = LTStransportsolver;
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, index)
            End = 0;
            
            % This is the outer loop
            obj.Converged = 0;
            obj.NLiter = 0;
            obj.NLiterLTS = 0;
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
                while obj.Converged == 0 && obj.itCount <= obj.MaxIter
                    % copy state to check outer convergence
                    State_old = status();
                    State_old.CopyProperties(ProductionSystem.Reservoir.State);
                    disp(newline);
                    disp(['Outer iteration: ', num2str(obj.itCount)]);
                    
                    Formulation.MatrixAssembler.ResetActiveInterfaces(DiscretizationModel);

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
%                     conservative = Formulation.CheckMassConservation(DiscretizationModel.ReservoirGrid);
%                     if ~conservative
%                         error('DARSim2 error: mass balance not respected');
%                     end
                    obj.BalanceTimer(obj.itCount) = toc(tstart2);
                    
                    % Save state before to compute the saturation
                    
                    State_iniTransp = status();
                    State_iniTransp.CopyProperties(ProductionSystem.Reservoir.State);

                    disp('Transport Solver');
                    disp('...............................................');
                    tstart3 = tic;
                    % 3.2 Solve Transport tentative global step
                    if obj.itCount == 1
                        obj.TransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
                    end
                    
                    obj.TransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                    obj.NLiter = obj.NLiter + obj.TransportSolver.itCount - 1;
                    obj.CFLGlobal = Formulation.ComputeCFLNumberTransport(DiscretizationModel, ProductionSystem, dt);
                    obj.TransportTimer(obj.itCount) = toc(tstart3);
                    if obj.TransportSolver.Converged == 0
                        disp('Transport solver failed to converge!');
                        obj.itCount = obj.MaxIter+1;
                    else
                        % 3.3 Compute active cells
                        obj.RefCellsSelector.SelectRefCells(ProductionSystem, DiscretizationModel.ReservoirGrid, Formulation);
                        
                        %if there is at least one active cell we need to
                        %recompute the new smaller system
                        if sum(sum(obj.RefCellsSelector.ActCells)) > 0
                            idxSummary = 1;
                            obj.ActCellsSummary = obj.RefCellsSelector.ActCells(:);
                            StateSum = status(); 
                            StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                            obj.StatesSummary = StateSum;
                            
                            DiscretizationModel.ReservoirGrid.ActiveTime = obj.RefCellsSelector.ActCells;
                            State_global = status();
                            State_global.CopyProperties(ProductionSystem.Reservoir.State);
                            dtRef = dt / obj.NofRef;
                            %Set initial values for the saturation
                            ProductionSystem.Reservoir.State.CopyProperties(State_iniTransp);
                           
                            % Compute the numerical fluxes used as boundary
                            % values between the accepted and rejected area.
                            obj.RefCellsSelector.SetActiveInterfaces(Formulation.MatrixAssembler, DiscretizationModel.ReservoirGrid)
                            obj.LTSTransportSolver.SystemBuilder.LTSBCEnforcer.SetCorrectActiveCells(obj.RefCellsSelector)
                            obj.LTSTransportSolver.SystemBuilder.LTSBCEnforcer.ComputeBoundaryValues(DiscretizationModel, Formulation, obj.RefCellsSelector);
                            
                          
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
                                ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                                
                                % save summary data 
                                obj.ActCellsSummary(:,idxSummary) = obj.RefCellsSelector.ActCells(:);
                                StateSum = status(); 
                                StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                                obj.StatesSummary(idxSummary) = StateSum;
                                idxSummary = idxSummary + 1;       
                                
                                obj.NLiterLTS = obj.NLiterLTS + obj.LTSTransportSolver.itCount - 1;
                                
                                if obj.LTSTransportSolver.Converged == 0
                                    disp('Sub ref not converged')
                                    obj.itCount = obj.MaxIter+1;
                                    break;
                                end
                                
                                ProductionSystem.SavePreviousState();
                                %for each subref
                                obj.LTSTransportTimer(obj.itCount) = obj.LTSTransportTimer(obj.itCount) + toc(tstart4);
                                disp('...............................................');
                            end
                            
                            disp('...............................................');
                            obj.CFLLocal = obj.LTSTransportSolver.SystemBuilder.LTSBCEnforcer.ComputeCFLNumberLTS(DiscretizationModel, ProductionSystem, dtRef, Formulation);

                            obj.LTSTransportSolver.SynchronizeProperties(ProductionSystem, State_global, obj.RefCellsSelector);
                            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                        else
                            obj.LTSTransportSolver.Converged = 1;
                            obj.LTSTransportTimer(obj.itCount) = 0;
                            obj.CFLLocal = 0;
                            obj.NLiterLTS = 0;
                            DiscretizationModel.ReservoirGrid.ActiveTime = zeros(DiscretizationModel.ReservoirGrid.N,1);
                        end
                        if obj.LTSTransportSolver.Converged ~= 0
                            if obj.MaxIter == 1 
                                % if there's only 1 outer-loop this check does not make sense
                                obj.Converged = 1;
                                obj.itCount = obj.itCount + 1;
                            else
                                obj.Converged = obj.ConvergenceChecker.Check(ProductionSystem.Reservoir.State, State_old);
                                obj.itCount = obj.itCount + 1;
                            end
                        else
                            obj.Converged = 0;
                        end
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
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
            obj.TimeStepSelector.UpdateSequential(dt, obj.itCount - 1, obj.TransportSolver.itCount - 1, obj.Chops);
            
        end
        function UpdateSummary(obj, Summary, Wells, Ndt, dt)
            Summary.CouplingStats.SaveStats(Ndt, obj.itCount - 1, obj.NLiter, obj.CFLGlobal, obj.NLiterLTS, obj.CFLLocal);
            Summary.CouplingStats.SaveTimers(Ndt, obj.PressureTimer, obj.BalanceTimer, obj.TransportTimer);
            Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
        end
    end
end