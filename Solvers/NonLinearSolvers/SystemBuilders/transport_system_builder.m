%  Transport System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%Last modified: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef transport_system_builder < system_builder
    properties
    end
    methods 
        function Residual = BuildResidual(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt)
           Residual = Formulation.BuildTransportResidual(ProductionSystem, FluidModel, DiscretizationModel, dt);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, FluidModel, DiscretizationModel, dt)
            Jacobian = Formulation.BuildTransportMatrix(ProductionSystem, FluidModel, DiscretizationModel, dt);
        end
    end
end