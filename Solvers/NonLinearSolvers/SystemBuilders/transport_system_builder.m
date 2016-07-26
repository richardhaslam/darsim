%  Transport System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef transport_system_builder < system_builder
    properties
    end
    methods
        function ComputePropertiesAndDerivatives(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel)
            Formulation.UpdateFractionalFlow(ProductionSystem, FluidModel);
            Formulation.dfdS(ProductionSystem, FluidModel);
        end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt)
           Residual = Formulation.BuildTransportResidual(ProductionSystem, DiscretizationModel, dt, obj.State);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt)
            Jacobian = Formulation.BuildTransportJacobian(ProductionSystem, DiscretizationModel, dt);
        end
        function delta = UpdateState(obj, delta, ProductionSystem, Formulation, FluidModel)
            Formulation.UpdateSaturation(ProductionSystem.Reservoir.State, delta, FluidModel);
        end
    end
end