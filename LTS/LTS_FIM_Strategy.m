% LTS FIM strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_FIM_Strategy < FIM_Strategy
    properties
        RefCellsSelector
    end
    methods
        function obj = LTS_FIM_Strategy(name, NONLinearSolver)
            obj@FIM_Strategy(name, NONLinearSolver);
            obj.RefCellsSelector = RefCellsSelector_BL();
        end
        function [dt, End] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation)
            % Initialise
            dt = obj.TimeStepSelector.ChooseTimeStep();
            dtf = dt / 10;
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
            
            %% 2. Compute flux b.c. for dtf zone
            Refinement = obj.ChooseRefZone(); % select indexes of ref zone cells
            
            %% 3. Solve fine time-step zone
            % Reset the state of refine zones to time-step n to solve again
            t = 0;
            while t < dt
                Formulation.Reset(); % It's only important for compositional
                % 3.A Set Up non-linear solver
                obj.NLSolver.SetUpLinearSolver(ProductionSystem, DiscretizationModel);
                obj.NLSolver.SetUp(Formulation, ProductionSystem, FluidModel, DiscretizationModel, dtf);
                % 3.B NL solver call
                obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dtf);
                % 3.C Update local time
                t = t + dtf;
            end
        end
        function Cells = ChooseRefZone(obj, ProductionSystem)
            Cells = obj.RefCellsSelector(ProductionSystem);
        end
        function ComputeBcRefZone(obj, Refinement)
            % If the neighbour of a fine cell is in the coarse zone I need to
            % compute the flux at n+1 and use it as a Neumann b.c.
            Cells = find(Refinement);
            for c = Cells
                for nc = Grid.Neighbours(c, :)
                    if Ref(nc)
                        % Flux for b.c.
                        obj.LTSBcEvaluator()
                    end
                end
            end
            
        end
    end
end

