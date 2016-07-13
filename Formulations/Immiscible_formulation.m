% Immiscible Formulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Immiscible_formulation < formulation
    properties
        
    end
    methods
        function BuildResidual(obj)
        end
        function BuildJacobian(obj)
        end
        function UpdateStatus(obj, ProductionSystem, delta)
        end
        function BuildPressureSystem(obj)
        end
        function BuildTransportSystem(obj)
        end
    end
end