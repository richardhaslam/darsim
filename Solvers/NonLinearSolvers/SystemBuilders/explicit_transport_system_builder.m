%  Explicit transport system Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 26 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef explicit_transport_system_builder < system_builder
    properties
    end
    methods
        function ComputePropertiesAndDerivatives(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel)
            Formulation.UpdateFractionalFlow(ProductionSystem, FluidModel);
        end
        function Residual = BuildResidual(obj)
            % not used: it's a virtual call
        end
        function Jacobian = BuildJacobian(obj)
            % not used: it's a virtual call
        end
    end
end