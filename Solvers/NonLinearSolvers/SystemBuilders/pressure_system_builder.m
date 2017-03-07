%  Pressure System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%Last modified: 4 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef pressure_system_builder < system_builder
    properties
    end
    methods 
        function ComputePropertiesAndDerivatives(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel)
            Formulation.ComputePropertiesAndDerivatives(ProductionSystem, FluidModel);
            Formulation.UpWindAndPhaseRockFluxes(DiscretizationModel.ReservoirGrid, FluidModel.Phases, ProductionSystem.Reservoir.State);
        end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt)
           Residual = Formulation.BuildPressureResidual(ProductionSystem, DiscretizationModel, dt, obj.State);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt)
            Jacobian = Formulation.BuildPressureMatrix(ProductionSystem, DiscretizationModel, dt);
        end
        function UpdateState(obj, delta, ProductionSystem, Formulation, FluidModel, DiscretizationModel)
            % Update Reservoir State
            Formulation.UpdatePressure(delta, ProductionSystem, FluidModel, DiscretizationModel);
            % UpdateWells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end
    end
end