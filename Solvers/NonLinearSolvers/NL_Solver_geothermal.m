% NL solver class for geothermal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NL_Solver_geothermal < NL_Solver
properties
end
methods
    function SetUp(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt)
        % Save initial state
        obj.SystemBuilder.SaveInitialState(ProductionSystem, Formulation);
        
        % Use CPR for Pressure (only for geothermal)
        switch(FluidModel.name)
            case {'Geothermal_1T','Geothermal_2T'}
                obj.SystemBuilder.ComputePropertiesAndDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
                Formulation.ConstrainedPressureResidual(FluidModel, ProductionSystem, DiscretizationModel, dt, obj.SystemBuilder.State);
                Formulation.ConstrainedTemperatureResidual(FluidModel, ProductionSystem, DiscretizationModel, dt, obj.SystemBuilder.State);
                obj.SystemBuilder.ComputePropertiesAndDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
        end
    end
    function SetUpLinearSolver(obj, ProductionSystem, DiscretizationModel)
        % Set up the linear solver
        obj.LinearSolver.SetUp(ProductionSystem, DiscretizationModel, 0); % We may need the residual for ADM but for now let's pass a 0.
    end
    function SetUpLinearSolverCoarse(obj, ProductionSystem, DiscretizationModel)
        obj.LinearSolver.LTS_SetUpCoarse(ProductionSystem, DiscretizationModel, obj.Residual);
    end
end
end