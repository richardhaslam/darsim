%  Natural variable Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NaturalVar_formulation < fim_formulation
    properties
        
    end
    methods
        function BuildResidual(obj)
        end
        function BuildJacobian(obj)
        end
        function UpdateState(obj, delta, Status, FluidModel)
            
        end
        function BuildPressureSystem(obj)
        end
        function BuildTransportSystem(obj)
        end
    end
end