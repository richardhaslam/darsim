% LTS FIM strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_FIM_Strategy < FIM_Strategy
    properties
        RefCellsSelector
        LTSNLSolver
        LTSNLSolverTimer
        NofRef = 5;
        LTSComplexity
    end
    methods
        function obj = LTS_FIM_Strategy(name, NONLinearSolver, LTSNONLSolver)
            obj@FIM_Strategy(name, NONLinearSolver);
            obj.LTSNLSolver = LTSNONLSolver;
            obj.RefCellsSelector = RefCellsSelector();
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
            % Initialise
            dt = obj.TimeStepSelector.ChooseTimeStep();
            dtRef = dt / obj.NofRef;
            obj.Converged = 0;
            obj.chops = 0;
            obj.LTSComplexity = 0;
            End = 0;
            
            %% 1. Solve with coarse time-step (predictor)
            % 1.1 Set Up non-linear solverLTSNONLSolver
            ProductionSystem.SavePreviousState(); % Save state of current time-step
            Formulation.MatrixAssembler.ResetActiveInterfaces(DiscretizationModel);
            obj.NLSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
            obj.NLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
            
            % 1.2 NL solver call
            disp(newline);
            disp('Global dt step');
            disp('...............................................');
            obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            
            %             figure(1)
            %             plot(ProductionSystem.Reservoir.State.Properties('S_1').Value, 'red');
            %             hold on
            %             figure(2)
            %             plot(ProductionSystem.Reservoir.State.Properties('P_1').Value, 'red')
            %             hold on
            
            if obj.NLSolver.Converged == 0
                disp('Global step did not converge')
                obj.Converged = 0;
                return
            else
                obj.Converged = 1;
            end
            
            
            if obj.TimeStepSelector.Index > 5
                %% 2 Select refined cells and compute fluxes
                % 2.0 Select refined cells
                Formulation.ComputeTotalFluxes(ProductionSystem, DiscretizationModel);
                % 2.1 Select refined cells
                obj.RefCellsSelector.SelectRefCells(ProductionSystem, DiscretizationModel.ReservoirGrid, Formulation);
                obj.RefCellsSelector.SetActiveInterfaces(Formulation.MatrixAssembler, DiscretizationModel.ReservoirGrid);
                DiscretizationModel.ReservoirGrid.ActiveTime = obj.RefCellsSelector.ActCells;
                % 2.2 Compute the numerical fluxes used as boundary values between the coarse and fine timesteps areas.
                obj.LTSNLSolver.SystemBuilder.LTSBCEnforcer.ComputeBoundaryFluxes(DiscretizationModel, Formulation, obj.RefCellsSelector.ActCells);
                
                %                 figure(3)
                %                 plot(obj.RefCellsSelector.ActCells, 'green')
                %
                
                %% 3. Solve fine time-step zone
                % 3.1 Reset the state of small dt zones to time-step n to solve again
                % I am not sure it is needed but probably yes
                obj.RefCellsSelector.ResetInitialState(ProductionSystem.Reservoir.State, ProductionSystem.Reservoir.State_old);
                
                % 3.2 Solve with small dt until sync is reached
                % 3.2.A Set Up non-linear solver
                obj.LTSNLSolver.LinearSolver.SetUp(ProductionSystem, DiscretizationModel,...
                    obj.RefCellsSelector.ActCells, FluidModel.NofPhases); % only once for now
                obj.LTSNLSolverTimer = 0;
                itSubRef = 1;
                
                disp(newline);
                obj.LTSComplexity = (obj.NLSolver.itCount-1) * length(obj.RefCellsSelector.ActCells);
                while itSubRef <= obj.NofRef && obj.Converged
                    disp(['SubRef step: ', num2str(itSubRef)]);
                    disp('...............................................');
                    tstart2 = tic;
                    
                    % 3.2.B Set Up non-linear solver
                    obj.LTSNLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dtRef);
                    % 3.2.C NL solver call
                    obj.LTSNLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dtRef);
                    
                    obj.LTSComplexity = obj.LTSComplexity + (obj.LTSNLSolver.itCount - 1)*sum(obj.RefCellsSelector.ActCells);
                    
                    if obj.LTSNLSolver.Converged == 0
                        disp(['Sub ref ', num2str(itSubRef),' did not converge']);
                        obj.Converged = 0;
                        break;
                    end
                    itSubRef = itSubRef + 1;
                    
                    %for each subref
                    obj.LTSNLSolverTimer = obj.LTSNLSolverTimer + toc(tstart2);
                    disp('...............................................');
                end
                
                %                 figure(1)
                %                 plot(ProductionSystem.Reservoir.State.Properties('S_1').Value, 'blue');
                %                 hold on
                %                 figure(2)
                %                 plot(ProductionSystem.Reservoir.State.Properties('P_1').Value, 'blue')
                %                 hold on
                %                 drawnow
            end
            obj.CFL = 0;
        end
        function Summary = UpdateSummary(obj, Summary, Wells, Ndt, dt)
            %% Stats, timers and Injection/Production data
            Summary.CouplingStats.SaveStats(Ndt, obj.NLSolver.itCount-1, obj.chops, obj.CFL);
            Summary.CouplingStats.SaveTimers(Ndt, obj.NLSolver.TimerConstruct, obj.NLSolver.TimerSolve, obj.NLSolver.TimerInner);
            Summary.CouplingStats.SaveLevelsComplexities(Ndt, obj.LTSComplexity);
            Summary.SaveWellsData(Ndt+1, Wells.Inj, Wells.Prod, dt);
        end
    end
end

