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
    MaxIter
    Converged
    itCount
    LinearSolver
    ConvergenceChecker
end
properties (Access = private)
    Residual
    Jacobian
    delta
end
methods
    function AddConvergenceChecker(obj, convcheck)
        obj.ConvergenceChecker = convcheck;
    end
    function [ProductionSystem, Summary] = Solve(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt, Summary)        
        % Initialise objects for new NL Solve
        obj.TimerConstruct = zeros(obj.MaxIter,1);
        obj.TimerSolve = zeros(obj.MaxIter, 1);
        obj.TimerInner = zeros(obj.MaxIter, 1);
        obj.itCount = 1;
        obj.Converged = 0;
        obj.CompConverged = 0;
        
        % Compute residual
        obj.BuildResidual(ProductionSystem, FluidModel, DiscretizationModel, dt);
        
        % Print some info
        disp(['Initial residual norm: ', num2str(norm(obj.Residual, inf))]);
        disp('');
        disp('        ||Residual||   ||delta p||   ||delta S||');
        
        % NEWTON-RAPHSON LOOP
        while ((obj.Converged==0)  && (obj.itCount <= obj.MaxIter))
            
            % 1. Build Jacobian Matrix for nu+1: everything is computed at nu
            start1 = tic;
            obj.BuildJacobian(ProductionSystem, FluidModel, DiscretizationModel, dt);
            obj.TimerConstruct(obj.itCount) = toc(start1);
            
            % 2. Solve full system at nu+1: J(nu)*Delta(nu+1) = -Residual(nu)
            start2 = tic;
            obj.SolveLinearSystem();
            obj.TimerSolve(obj.itCount) = toc(start2);
            
            % 3. Update Solution
            start3 = tic;
            ProductionSystem = Formulation.UpdateSolution(ProductionSystem);
            obj.TimerInner(obj.itCount) = toc(start3);
            
            % 5. Compute residual
            obj.BuildResidual(ProductionSystem, FluidModel, DiscretizationModel, dt);
            
            % 6. Check NonLinear convergence
            obj.CheckConvergence();
            
            obj.itCount = obj.itCount + 1;
        end
    end
    function BuildResidual(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel, dt)
        obj.Residual = Formulation.BuildResidual(ProductionSystem, FluidModel, DiscretizationModel, dt);
    end
    function SolveLinearSystem(obj)
        obj.delta = obj.LinearSolver.Solve(obj.Jacobian, obj.Residual);
    end
    function BuildJacobian(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel, dt)
        obj.Jacobian = Formulation.BuildJacobian(ProductionSystem, FluidModel, DiscretizationModel, dt);
    end
    function CheckConvergence(obj)
        obj.Converged = obj.ConvergenceChecker.Check();
    end
end
end