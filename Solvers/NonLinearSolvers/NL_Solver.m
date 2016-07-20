% NL solver base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 14 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NL_Solver < handle
properties
    SystemBuilder
    MaxIter
    Converged
    itCount
    LinearSolver
    ConvergenceChecker
    TimerConstruct
    TimerSolve
    TimerInner
end
properties (Access = private)
    Residual
    Jacobian
    Delta
end
methods
    function obj = NL_Solver()
        obj.Converged = 0;
    end
    function AddConvergenceChecker(obj, convcheck)
        obj.ConvergenceChecker = convcheck;
    end
    function ProductionSystem = Solve(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt)        
        % Initialise objects for new NL Solve
        obj.TimerConstruct = zeros(obj.MaxIter,1);
        obj.TimerSolve = zeros(obj.MaxIter, 1);
        obj.TimerInner = zeros(obj.MaxIter, 1);
        obj.itCount = 1;
        
        % Save initial State
        obj.SystemBuilder.SaveInitialState(ProductionSystem.Reservoir.State);
        
        % Update Derivatives
        Formulation = obj.SystemBuilder.ComputeDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
        % Compute residual
        obj.BuildResidual(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
        
        % Print some info
        disp(['Initial residual norm: ', num2str(norm(obj.Residual, inf))]);
        disp('');
        disp('        ||Residual||   ||delta p||   ||delta S||');
        
        % NEWTON-RAPHSON LOOP
        while ((obj.Converged==0)  && (obj.itCount <= obj.MaxIter))
            
            % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
            start1 = tic;
            obj.BuildJacobian(ProductionSystem, Formulation, DiscretizationModel, dt);
            obj.TimerConstruct(obj.itCount) = toc(start1);
            
            % 2. Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
            start2 = tic;
            obj.SolveLinearSystem();
            obj.TimerSolve(obj.itCount) = toc(start2);
            
            % 3. Update Solution
            start3 = tic;
            ProductionSystem.UpdateState(obj.Delta, Formulation, FluidModel);
            obj.TimerInner(obj.itCount) = toc(start3);
            
            % 4. Update Derivatives
            Formulation = obj.SystemBuilder.ComputeDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
            % 5. Compute residual
            obj.BuildResidual(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            
            % 6. Check NonLinear convergencel
            obj.CheckConvergence(DiscretizationModel, ProductionSystem);
            
            obj.itCount = obj.itCount + 1;
        end
    end
    function BuildResidual(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt)
        obj.Residual = obj.SystemBuilder.BuildResidual(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt); 
    end
    function SolveLinearSystem(obj)
        obj.Delta = obj.LinearSolver.Solve(obj.Jacobian, -obj.Residual);
    end
    function BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt)
        obj.Jacobian = obj.SystemBuilder.BuildJacobian(ProductionSystem, Formulation, DiscretizationModel, dt);
    end
    function CheckConvergence(obj, DiscretizationModel, ProductionSystem)
        obj.Converged = obj.ConvergenceChecker.Check(obj.itCount, obj.Residual, obj.Delta, DiscretizationModel, ProductionSystem.Reservoir.State.p);
    end
end
end