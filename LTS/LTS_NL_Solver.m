% NL solver base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_NL_Solver < NL_Solver
properties
    LTSActiveCells
end
methods
    function BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt)
        obj.Residual = obj.SystemBuilder.BuildResidual(ProductionSystem, DiscretizationModel, Formulation, dt);
    end
    function BuildJacobian(obj, ProductionSystem, DiscretizationModel, Formulation, dt)
        obj.Residual = obj.SystemBuilder.BuildResidual(ProductionSystem, DiscretizationModel, Formulation, dt);
    end
    function UpdateState(obj, ProductionSystem, Formulation, FluidModel, DiscretizationModel)
        obj.Delta = obj.SystemBuilder.UpdateState(obj.Delta, ProductionSystem, Formulation, FluidModel, DiscretizationModel);
    end
    function SetUp(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt)
        % 1. Save initial state
        obj.SystemBuilder.SaveInitialState(ProductionSystem, Formulation);
        
        % 2. Computes initial residual
        % Update Derivatives
        obj.SystemBuilder.ComputePropertiesAndDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
        % Compute residual
        obj.BuildResidual(ProductionSystem, DiscretizationModel, Formulation, dt);
    end
    function SetUpLinearSolver(obj, ProductionSystem, DiscretizationModel)
        % Set up the linear solver
        obj.LinearSolver.SetUp(ProductionSystem, DiscretizationModel, obj.Residual);
    end
end
end