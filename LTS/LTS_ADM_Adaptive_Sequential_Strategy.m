classdef LTS_ADM_Adaptive_Sequential_Strategy < LTS_Adaptive_Sequential_Strategy
    properties
        time_ref = 2;
        LTS_iters;
    end
    methods
        function obj = LTS_ADM_Adaptive_Sequential_Strategy(name)
            obj@ LTS_Adaptive_Sequential_Strategy(name);
            obj.RefCellsSelector = RefCellsSelector();
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
            End = 0;
            obj.LTS_iters = [];
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
                    %obj.PressureSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
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
                    %It is not correct
                    %                     conservative = Formulation.CheckMassConservation(DiscretizationModel.ReservoirGrid);
                    %                     if ~conservative
                    %                         error('DARSim2 error: mass balance not respected');
                    %                     end
                    obj.BalanceTimer(obj.itCount) = toc(tstart2);
                    
                    % Save state before to compute the saturation
                    State_iniTransp = status();
                    State_iniTransp.CopyProperties(ProductionSystem.Reservoir.State);
                    
                    obj.TransportSolver.SetUpLinearSolverCoarse(ProductionSystem, DiscretizationModel);
                    %% 3. Solve transport
                      % we do not se it in order to impose fix global steps
%                     if obj.itCount == 1 && obj.Chops == 0
%                         dt = obj.TimeStepSelector.StableTimeStep(ProductionSystem, DiscretizationModel, FluidModel, Formulation.Utot);
%                     end
                    disp('Transport Solver');
                    disp('...............................................');
                    tstart3 = tic;
                    % 3.2 Solve Transport tentative global step on the
                    % coarse grid
                    if obj.itCount == 1
                        obj.TransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
                    end
                    obj.TransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                    obj.NLiter = obj.NLiter + (obj.TransportSolver.itCount - 1)*DiscretizationModel.ADMGrid.Ntot;
                    obj.CFLGlobal = Formulation.ComputeCFLNumberTransport(DiscretizationModel, ProductionSystem, dt);
                    obj.TransportTimer(obj.itCount) = toc(tstart3);
                    if obj.TransportSolver.Converged == 0
                        disp('Transport solver failed to converge!');
                        obj.itCount = obj.MaxIter+1;
                    else
                        obj.TransportSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
                        obj.LTSTransportSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);% store the global state -> sync of cells at the
                        % end
                        State_global = status();
                        State_global.CopyProperties(ProductionSystem.Reservoir.State);
                        obj.StateGlobalVec = State_global;
                        
                        NofLevel = DiscretizationModel.maxLevel; % num of Level that we need to comp
                        lev = 1;
                        dtRef = dt;
                        % Compute active cells on the basis of the active cells on ADM grid
                        obj.RefCellsSelector.ComputeActiveCells(DiscretizationModel, lev);
                        
                        
                        % store it in the vector and compute BC latent
                        % cells
                        RefCells = RefCellsSelector();
                        RefCells.CopyCellsSelected(obj.RefCellsSelector);
                        obj.RefCellsSelectorVec = RefCells;
                        
                        obj.RefCellsSelectorVec.ComputeBoundaryValues(DiscretizationModel, Formulation);
                        
                        dtGlob = dtRef;
                        dtRef = dtRef / obj.time_ref;
                        sum_dtLoc = zeros(factorial(NofLevel) * obj.time_ref);
                        
                        idxSummary = 1;
                        % to store local time step
                        obj.ActCellsSummary = obj.RefCellsSelector.ActCells(:);
                        StateSum = status();
                        StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                        obj.StatesSummary = StateSum;
                        
                        %Set initial values for the saturation
                        ProductionSystem.Reservoir.State.CopyProperties(State_iniTransp);
                        ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                        
                        Newton_IniGuess = status();
                        Newton_IniGuess.CopyProperties(ProductionSystem.Reservoir.State);
                        temp = Newton_IniGuess.Properties('S_1');                        
                        temp2 = Newton_IniGuess.Properties('S_2');

                        S1_next = obj.StateGlobalVec(1).Properties('S_1').Value;
                        S2_next = obj.StateGlobalVec(1).Properties('S_2').Value;
                        S1_old = State_iniTransp.Properties('S_1').Value(:);
                        S2_old = State_iniTransp.Properties('S_2').Value(:);
                        temp.Value = S1_old + (dtRef)/dtGlob * (S1_next - S1_old);
                        temp2.Value = S2_old + (dtRef)/dtGlob * (S2_next - S2_old);

                        while sum_dtLoc(lev) < dtGlob(lev)
                            % we sum up all the timer for the refinemets
                            obj.LTSTransportTimer(obj.itCount) = 0;
                            disp('Transport Solver Sub-rebinements');
                            disp('...............................................');
                            
                            State_iniTransp = status();
                            State_iniTransp.CopyProperties(ProductionSystem.Reservoir.State);
                            
                            disp('...............................................');
                            disp(['SubRef step: ', num2str(lev)]);
                            disp('...............................................');
                            % Set up for LinearSolver R and P
                            if lev == NofLevel
                                obj.LTSTransportSolver.SetUpRP_LTS_ADM(DiscretizationModel, obj.RefCellsSelectorVec(lev).ActCells(:), lev)  
                                obj.LTSTransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dtRef, obj.RefCellsSelectorVec(lev));
                                 
                                obj.LTSTransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dtRef, obj.RefCellsSelectorVec(lev));
                                obj.LTS_iters = [obj.LTS_iters (obj.LTSTransportSolver.itCount-1)*obj.RefCellsSelectorVec(lev). NumberOfActiveCells(DiscretizationModel, lev)];
                            
                                obj.CFLLocal = Formulation.ComputeCFLNumberTransportLTS(DiscretizationModel, ProductionSystem, dtRef, obj.RefCellsSelectorVec(lev));
                                obj.NLiterLTS = obj.NLiterLTS + obj.LTSTransportSolver.itCount - 1;
                            else 
                                obj.LTSTransportSolver.SetUpRP_LTS_ADM(DiscretizationModel, obj.RefCellsSelectorVec(lev).ActCells(:), lev)
                                obj.LTSTransportSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dtRef, obj.RefCellsSelectorVec(lev));
                                
                                ProductionSystem.Reservoir.State.CopyProperties(Newton_IniGuess);
                                                                
                                obj.LTSTransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dtRef, obj.RefCellsSelectorVec(lev));
                                obj.LTS_iters = [obj.LTS_iters (obj.LTSTransportSolver.itCount-1)*obj.RefCellsSelectorVec(lev). NumberOfActiveCells(DiscretizationModel, lev)];
                            
                                obj.CFLLocal = Formulation.ComputeCFLNumberTransportLTS(DiscretizationModel, ProductionSystem, dtRef, obj.RefCellsSelectorVec(lev));
                                obj.NLiterLTS = obj.NLiterLTS + obj.LTSTransportSolver.itCount - 1;
                            end
                                                          
                            if lev == DiscretizationModel.maxLevel
                                sum_dtLoc(lev) = sum_dtLoc(lev) + dtRef;
                                
                                % save summary data
                                obj.ActCellsSummary(:,idxSummary) = obj.RefCellsSelectorVec(lev).ActCells(:);
                                StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                                obj.StatesSummary(idxSummary) = StateSum;
                                idxSummary = idxSummary + 1;
                                
                                
                                while abs(sum_dtLoc(lev) - dtGlob(lev)) < eps && lev > 1
                                    
                                    obj.LTSTransportSolver.SynchronizeProperties(ProductionSystem, obj.StateGlobalVec(lev), obj.RefCellsSelectorVec(lev));
                                    ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                                    lev = lev - 1;
                                    sum_dtLoc(lev) = sum_dtLoc(lev) + sum_dtLoc(lev + 1);
                                    dtRef = obj.time_ref * dtRef;
                                    
                                    idxSummary = idxSummary - 1;                                 % save summary data
                                    obj.ActCellsSummary(:,idxSummary) = obj.RefCellsSelectorVec(lev).ActCells(:);
                                    StateSum = status();
                                    StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                                    obj.StatesSummary(idxSummary) = StateSum;
                                    idxSummary = idxSummary + 1;
                                    temp = Newton_IniGuess.Properties('S_1');
                                    temp2 = Newton_IniGuess.Properties('S_2');
                                    S1_next = obj.StateGlobalVec(1).Properties('S_1').Value;
                                    S2_next = obj.StateGlobalVec(1).Properties('S_2').Value;
                                    S1_old = State_iniTransp.Properties('S_1').Value(:);
                                    S2_old = State_iniTransp.Properties('S_2').Value(:);
                                    if (sum_dtLoc(lev)+ dtRef - dtGlob(lev))< eps
                                        temp.Value = S1_old + (sum_dtLoc(lev)+ dtRef)/dtGlob(lev) * (S1_next - S1_old);
                                        temp2.Value = S2_old + (sum_dtLoc(lev)+ dtRef)/dtGlob(lev) * (S2_next - S2_old);
                                    else
                                        temp.Value = S1_next;
                                        temp2.Value = S2_next;
                                    end
                                end
                            else
                                obj.RefCellsSelector.ComputeActiveCells(DiscretizationModel, lev+1);
                                if lev == DiscretizationModel.maxLevel - 1
                                    obj.RefCellsSelector.ActCellCheckError(ProductionSystem, DiscretizationModel.ReservoirGrid, Formulation);
                                end
                                DiscretizationModel.ReservoirGrid.ActiveTime = obj.RefCellsSelector.ActCells;

                                if sum(obj.RefCellsSelector.ActCells) ~= 0
                                    lev = lev + 1;
                                    State_global = status();
                                    State_global.CopyProperties(ProductionSystem.Reservoir.State);
                                    obj.StateGlobalVec(lev) = State_global;
                                    
                                    RefCells = RefCellsSelector();
                                    RefCells.CopyCellsSelected(obj.RefCellsSelector);
                                    obj.RefCellsSelectorVec(lev) = RefCells;
                                    obj.RefCellsSelectorVec(lev).ComputeBoundaryValuesSubRef(DiscretizationModel, Formulation, obj.RefCellsSelectorVec(lev-1));
                                    
                                    %Set initial values for the saturation
                                    ProductionSystem.Reservoir.State.CopyProperties(State_iniTransp);
                                    ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                                    
                                    dtGlob(lev) = dtRef;
                                    dtRef = dtRef / obj.time_ref;
                                    sum_dtLoc(lev) = 0;
                                else
                                    sum_dtLoc(lev) = sum_dtLoc(lev) + dtRef;
                                    
                                    % save summary data
                                    obj.ActCellsSummary(:,idxSummary) = obj.RefCellsSelectorVec(lev).ActCells(:);
                                    StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                                    obj.StatesSummary(idxSummary) = StateSum;
                                    idxSummary = idxSummary + 1;
                                end
                            end
                        end
                        
                        disp('...............................................');
                        ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
                        obj.LTSTransportSolver.SynchronizeProperties(ProductionSystem, obj.StateGlobalVec(lev), obj.RefCellsSelectorVec(lev));

                        idxSummary = idxSummary - 1;                                 % save summary data
                        obj.ActCellsSummary(:,idxSummary) = ones(size(obj.RefCellsSelectorVec(lev).ActCells(:)));
                        StateSum = status();
                        StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                        obj.StatesSummary(idxSummary) = StateSum;
                        
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
