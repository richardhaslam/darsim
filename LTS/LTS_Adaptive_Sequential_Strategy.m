% LTS Sequential strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author:Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_Adaptive_Sequential_Strategy < LTS_Sequential_Strategy
    properties
        RefCellsSelectorVec
        StateGlobalVec
    end
    methods
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
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
            %dt = 5 * 60 * 60 * 24;
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
                while obj.Converged == 0 && obj.itCount <= obj.MaxIter
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
      
                    %% 3. Solve transport
                    if obj.itCount == 1 && obj.Chops == 0
                        dt = obj.TimeStepSelector.StableTimeStep(ProductionSystem, DiscretizationModel, FluidModel, Formulation.Utot);
                    end
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
                        
                        % if there is at least one active cell we need to
                        % recompute the new smaller system
                        if sum(sum(obj.RefCellsSelector.ActCells)) > 0 && obj.CFLGlobal > 1
                            % useful to store informations for the sub-ref
                            % steps; 
                            idxSummary = 1;
                            obj.ActCellsSummary = obj.RefCellsSelector.ActCells(:);
                            StateSum = status(); 
                            StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                            obj.StatesSummary = StateSum;
                            
                            % to count the refinement sub-steps;
                            itRef = 1;
                            % vector contains refCellst for
                            % each level of sub-refinement
                            RefCells =  RefCellsSelector();
                            RefCells.CopyCellsSelected(obj.RefCellsSelector)
                            obj.RefCellsSelectorVec = RefCells;
                            % a the moment we save just the active comp of the
                            % first sub-ref 
                            DiscretizationModel.ReservoirGrid.ActiveTime = obj.RefCellsSelector.ActCells;
                            State_global = status(); 
                            State_global.CopyProperties(ProductionSystem.Reservoir.State);
                            obj.StateGlobalVec = State_global;
                            
                            dtRef = dt/2;
                            dtGlob = dt;
                            sum_dtLoc(itRef) = 0;                           
                            
                            %Set initial values for the saturation
                            ProductionSystem.Reservoir.State.CopyProperties(State_iniTransp);
                            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                            
                            % Compute the numerical fluxes used as boundary
                            % values between the accepted and rejected area.
                            obj.RefCellsSelectorVec(itRef).ComputeBoundaryValues(DiscretizationModel, Formulation);
                            
                            % we sum up all the time for the refinemets
                            obj.LTSTransportTimer(obj.itCount) = 0;
                            disp('Transport Solver Sub-rebinements');
                           
                            while sum_dtLoc(itRef) < dtGlob(itRef)
                                % Solve the transport eq for the subref cells
                                disp('...............................................');
                                disp(['SubRef step: ', num2str(itRef)]);
                                disp('...............................................');
                                
                                State_iniTransp = status();
                                State_iniTransp.CopyProperties(ProductionSystem.Reservoir.State);
                                
                                obj.LTSTransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dtRef, obj.RefCellsSelectorVec(itRef));
                                obj.LTSTransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dtRef, obj.RefCellsSelectorVec(itRef));
                                
                                % at the moment I am just storing the
                                % smallest one
                                obj.CFLLocal = Formulation.ComputeCFLNumberTransportLTS(DiscretizationModel, ProductionSystem, dtRef, obj.RefCellsSelectorVec(itRef));

                                obj.NLiterLTS = obj.NLiterLTS + obj.LTSTransportSolver.itCount - 1;
                                                                
                                if obj.LTSTransportSolver.Converged == 0
                                    disp('Sub ref not converged')
                                    obj.itCount = obj.MaxIter+1;
                                    break;
                                end 
                                % check if we need a new sub-ref only if
                                % the CFL is > 1
                                obj.RefCellsSelector.SelectRefCells(ProductionSystem, DiscretizationModel.ReservoirGrid, Formulation);
                                
                                if sum(sum(obj.RefCellsSelector.ActCells)) > 0 && obj.CFLLocal > 1
                                    itRef = itRef + 1;
                                    RefCells =  RefCellsSelector();
                                    RefCells.CopyCellsSelected(obj.RefCellsSelector)
                                    obj.RefCellsSelectorVec(itRef) = RefCells;
                                    %Update the new boundary values
                                    obj.RefCellsSelectorVec(itRef).ComputeBoundaryValuesSubRef(DiscretizationModel, Formulation, obj.RefCellsSelectorVec(itRef-1)); 
                                    
                                    State_global = status(); 
                                    State_global.CopyProperties(ProductionSystem.Reservoir.State);
                                    obj.StateGlobalVec(itRef) = State_global;                            
                                    
                                    %Set initial values for the saturation
                                    ProductionSystem.Reservoir.State.CopyProperties(State_iniTransp);
                                    ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                                    
                                    dtGlob(itRef) = dtRef;
                                    dtRef = dtRef / 2;
                                    sum_dtLoc(itRef) = 0;
 
                                else
                                    % if we do not need a more sub-ref we can performe a new time step     
                                    sum_dtLoc(itRef) = sum_dtLoc(itRef) + dtRef; 
                                    
                                    % save summary data
                                    obj.ActCellsSummary(:,idxSummary) = obj.RefCellsSelectorVec(itRef).ActCells(:);
                                    StateSum = status(); 
                                    StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                                    obj.StatesSummary(idxSummary) = StateSum;
                                    idxSummary = idxSummary + 1;
                                    
                                    % check if we need to exit from a
                                    % sub-ref
                                    while abs(sum_dtLoc(itRef) - dtGlob(itRef)) < eps && itRef>1
                                        
                                        obj.LTSTransportSolver.SynchronizeProperties(ProductionSystem, obj.StateGlobalVec(itRef), obj.RefCellsSelectorVec(itRef));    
                                        ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                                        itRef = itRef - 1;
                                        sum_dtLoc(itRef) = sum_dtLoc(itRef) + sum_dtLoc(itRef + 1);
                                        dtRef = 2 * dtRef;
                                    end
                                end
                            end
                            
                            disp('...............................................');
                            obj.LTSTransportSolver.SynchronizeProperties(ProductionSystem, obj.StateGlobalVec(itRef), obj.RefCellsSelectorVec(itRef));
                            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
 
                        else
                            obj.LTSTransportSolver.Converged = 1;
                            obj.LTSTransportTimer(obj.itCount) = 0;
                            obj.CFLLocal = 0;
                            obj.NLiterLTS = 0;
                            obj.ActCellsSummary = 0;
                            obj.StatesSummary = 0;
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