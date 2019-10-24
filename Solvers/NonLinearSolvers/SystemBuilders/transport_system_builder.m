%  Transport System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef transport_system_builder < system_builder
    properties
    end
    methods
        function ComputePropertiesAndDerivatives(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel)
            Formulation.UpdateFractionalFlow(ProductionSystem, FluidModel);
            Formulation.UpdateCapillaryPressure(ProductionSystem, FluidModel);
            Formulation.ComputeDfDS(ProductionSystem, FluidModel);
        end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt)
           Residual = Formulation.BuildTransportResidual(ProductionSystem, DiscretizationModel, dt, obj.State);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt)
            Jacobian = Formulation.BuildTransportJacobian(ProductionSystem, DiscretizationModel, dt);
        end
        function delta = UpdateState(obj, delta, ProductionSystem, Formulation, FluidModel, DiscretizationModel, iter)
            delta = Formulation.UpdateSaturation(ProductionSystem, delta, FluidModel, DiscretizationModel, iter);
            % UpdateWells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end
        function SetUpSolutionChopper(obj, SolutionChopper, Formulation, ProductionSystem, N)
            x = Formulation.GetPrimaryPressure(ProductionSystem, N);
            SolutionChopper.DefineMaxDelta(x);
        end
    end
end