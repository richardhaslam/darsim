                                                                                                                                                                     % NL solver base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_NL_Solver < handle
properties
    SystemBuilder
    SolutionChopper
    MaxIter
    Converged
    itCount
    LinearSolver
    ConvergenceChecker
    TimerConstruct
    TimerSolve
    TimerInner
    Delta
    InitialTimeStep = true;
end
properties (Access = private)
    Residual
    Jacobian
%     InitialTimeStep = true;
end
methods
    function obj = LTS_NL_Solver()
        obj.Converged = 0;
        obj.SolutionChopper = solution_chopper();
    end
    function AddConvergenceChecker(obj, convcheck)
        obj.ConvergenceChecker = convcheck;
    end
    function Solve(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt, CellsSelected)
        % Initialise objects for new NL Solve
        obj.Converged = 0;
        obj.TimerConstruct = zeros(obj.MaxIter,1);
        obj.TimerSolve = zeros(obj.MaxIter, 1);
        obj.TimerInner = zeros(obj.MaxIter, 1);
        obj.itCount = 1;
        
        % 0. Computes initial residual
        % Update Derivatives
        obj.SystemBuilder.ComputePropertiesAndDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
        % Compute residual
        obj.BuildResidual(ProductionSystem, DiscretizationModel, Formulation, dt, CellsSelected);
        % Compute the first residual norm
        obj.ConvergenceChecker.ComputeFirstResidualNorm(obj.Residual, DiscretizationModel, obj.LinearSolver);
        
        % Print some info
        obj.ConvergenceChecker.PrintTitles();
   
        % NEWTON-RAPHSON LOOP
        while ((obj.Converged==0)  && (obj.itCount <= obj.MaxIter))
            
            % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
            start1 = tic;
            obj.BuildJacobian(ProductionSystem, Formulation, DiscretizationModel, dt, CellsSelected);
            obj.TimerConstruct(obj.itCount) = toc(start1);
            
            % 2. Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
            obj.SystemBuilder.SetUpSolutionChopper(obj.SolutionChopper, Formulation, ProductionSystem, DiscretizationModel);
            start2 = tic;
            obj.SolveLinearSystem(CellsSelected.ActCells);
            obj.TimerSolve(obj.itCount) = toc(start2);
            obj.Delta = obj.SolutionChopper.Chop(obj.Delta);
            
            % 3. Update Solution
            start3 = tic;
            obj.UpdateState(ProductionSystem, Formulation, FluidModel, DiscretizationModel);
            obj.TimerInner(obj.itCount) = toc(start3);
            
            % 4. Update Derivatives
            obj.SystemBuilder.ComputePropertiesAndDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
            
            % 5. Compute residual
            obj.BuildResidual(ProductionSystem, DiscretizationModel, Formulation, dt, CellsSelected);
            
            % 6. Check NonLinear convergencel
            obj.CheckConvergence(Formulation, DiscretizationModel, ProductionSystem);
            
            if obj.Converged == -1
                obj.Converged = 0;
                return
            end
            obj.itCount = obj.itCount + 1;
        end
    end
    function BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt, CellsSelected)
        obj.Residual = obj.SystemBuilder.BuildResidual(ProductionSystem, DiscretizationModel, Formulation, dt, CellsSelected);
    end
    function SolveLinearSystem(obj, ActCells)
        obj.Delta = obj.LinearSolver.Solve(obj.Jacobian, -obj.Residual, ActCells);
    end
    function BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt, CellsSelected)
        obj.Jacobian = obj.SystemBuilder.BuildJacobian(ProductionSystem, Formulation, DiscretizationModel, dt, CellsSelected);
    end
    function CheckConvergence(obj, Formulation, DiscretizationModel, ProductionSystem)
        obj.Converged = obj.ConvergenceChecker.Check(obj.itCount, obj.Residual, obj.Delta, Formulation, DiscretizationModel, ProductionSystem.Reservoir.State, obj.LinearSolver);
    end
    function UpdateState(obj, ProductionSystem, Formulation, FluidModel, DiscretizationModel)
        obj.Delta = obj.SystemBuilder.UpdateState(obj.Delta, ProductionSystem, Formulation, FluidModel, DiscretizationModel);
    end
    function SetUp(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt, CellsSelected)
        % 1. Save initial state
        obj.SystemBuilder.SaveInitialState(ProductionSystem, Formulation);
               
    end
    function SetUpLinearSolver(obj, ProductionSystem, DiscretizationModel)
        % Set up the linear solver
        obj.LinearSolver.SetUp(ProductionSystem, DiscretizationModel, obj.Residual);
    end
    
    function SynchronizeProperties(obj, ProductionSystem, State_global, CellsSelected)
        obj.SystemBuilder.SynchronizeProperties(ProductionSystem, State_global, CellsSelected);
    end
end
end