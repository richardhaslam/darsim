% Mass Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Mass_formulation < handle
    properties
        
    end
    methods
        function BuildResidual(obj)
        end
        function BuildJacobian(obj)
        end
        function UpdateStatus(obj, Status, delta)
        end
        function BuildPressureSystem(obj)
        end
        function BuildTransportSystem(obj)
        end
    end
end