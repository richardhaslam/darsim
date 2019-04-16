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
        NofRef = 10;
        NLiterLTS
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
            obj.NLiterLTS = 0;
            End = 0;
            
            %% 1. Solve with coarse time-step (predictor)
            % 1.1 Set Up non-linear solverLTSNONLSolver
            ProductionSystem.SavePreviousState(); % Save state of current time-step
            Formulation.MatrixAssembler.ResetAcitveInterfaces(DiscretizationModel);
            obj.NLSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
            obj.NLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
            
            % 1.2 NL solver call
            obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            
            %% 2 Select refined cells and compute fluxes
            % 2.0 Select refined cells
            Formulation.ComputeTotalFluxes(ProductionSystem, DiscretizationModel);
            % 2.1 Select refined cells
            obj.RefCellsSelector.SelectRefCells(ProductionSystem, DiscretizationModel.ReservoirGrid, Formulation);
            obj.RefCellsSelector.SetActiveInterfaces(Formulation.MatrixAssembler, DiscretizationModel.ReservoirGrid);
            % 2.2 Compute the numerical fluxes used as boundary values between the coarse and fine timesteps areas.
            obj.LTSNLSolver.SystemBuilder.LTSBCEnforcer.ComputeBoundaryFluxes(DiscretizationModel, Formulation, obj.RefCellsSelector.ActCells);
            
            %% 3. Solve fine time-step zone
            % 3.1 Reset the state of small dt zones to time-step n to solve again
            % I am not sure it is needed but probably yes
            obj.RefCellsSelector.ResetInitialState(ProductionSystem.Reservoir.State, ProductionSystem.Reservoir.State_old);
            
            % 3.2 Solve with small dt until sync is reached
            % 3.2.A Set Up non-linear solver
            obj.LTSNLSolver.LinearSolver.SetUp(ProductionSystem, DiscretizationModel,...
                obj.RefCellsSelector.ActCells, FluidModel.NofPhases); % only once for now
            obj.LTSNLSolverTimer = 0;
            for itSubRef = 1:obj.NofRef
                disp(['SubRef step: ', num2str(itSubRef)]);
                disp('...............................................');
                tstart2 = tic;
              
                % 3.2.B Set Up non-linear solver
                obj.LTSNLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dtRef);
                % 3.2.C NL solver call
                obj.LTSNLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dtRef);
                                
                obj.NLiterLTS = obj.NLiterLTS + obj.LTSNLSolver.itCount - 1;
                
                if obj.LTSNLSolver.Converged == 0
                    disp('Sub ref not converged')
                    obj.itCount = obj.MaxIter+1;
                    break;
                end
                
                %for each subref
                obj.LTSNLSolverTimer = obj.LTSNLSolverTimer + toc(tstart2);
                disp('...............................................');
            end
        end
    end
end

