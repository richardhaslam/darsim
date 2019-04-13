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
        NofRef = 10;
    end
    methods
        function obj = LTS_FIM_Strategy(name, NONLinearSolver)
            obj@FIM_Strategy(name, NONLinearSolver);
            obj.RefCellsSelector = RefCellsSelector();
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
            % Initialise
            dt = obj.TimeStepSelector.ChooseTimeStep();
            dtRef = dt / obj.NofRef;
            obj.Converged = 0;
            obj.chops = 0;
            End = 0;
            
            %% 1. Solve with coarse time-step (predictor)
            % 1.A Set Up non-linear solver
            obj.NLSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
            obj.NLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt);
            ProductionSystem.SavePreviousState(); % Save state of current time-step
            % 1.B NL solver call
            obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            
            %% 2 Select refined cells and compute fluxes
            % 2.0 Select refined cells
            Formulation.ComputeTotalFluxes(ProductionSystem, DiscretizationModel);
            % 2.1 Select refined cells
            obj.RefCellsSelector.SelectRefCells(ProductionSystem, DiscretizationModel.ReservoirGrid, Formulation);
            % 2.2 Compute the numerical fluxes used as boundary values between the coarse and fine timesteps areas.
            obj.RefCellsSelector.ComputeBoundaryValues(DiscretizationModel, Formulation);
                        
            %% 3. Solve fine time-step zone
            % 3.1 Reset the state of small dt zones to time-step n to solve again
            % I am not sure it is needed but probably yes
            
            % 3.2 Solve with small dt until sync is reached
            for i=1:obj.NofRef
                disp(['SubRef step: ', num2str(itSub)]);
                disp('...............................................');
                tstart2 = tic;
                
                Formulation.Reset(); % It's only important for compositional
                % 3.A Set Up non-linear solver
                obj.LTSNLSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
                obj.LTSNLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dtRef);
                % 3.B NL solver call
                obj.LTSNLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dtRef);
                
                % save summary data
                % obj.ActCellsSummary(:,idxSummary) = obj.RefCellsSelector.ActCells(:);
                StateSum = status();
                StateSum.CopyProperties(ProductionSystem.Reservoir.State);
                obj.StatesSummary(idxSummary) = StateSum;
                idxSummary = idxSummary + 1;
                
                obj.NLiterLTS = obj.NLiterLTS + obj.LTSNLSolver.itCount - 1;
                
                if obj.LTSNLSolver.Converged == 0
                    disp('Sub ref not converged')
                    obj.itCount = obj.MaxIter+1;
                    break;
                end
                
                ProductionSystem.SavePreviousState();
                %for each subref
                obj.LTSTimer(obj.itCount) = obj.LTSTransportTimer(obj.itCount) + toc(tstart2);
                disp('...............................................');
            end
        end
        function Cells = ChooseRefZone(obj, ProductionSystem)
            Cells = obj.RefCellsSelector(ProductionSystem);
        end
        
    end
end

