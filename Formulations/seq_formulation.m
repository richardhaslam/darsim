%  Sequential Formulation base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%Last modified: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef seq_formulation < handle
    properties
        NumOfEquations
    end
    methods 
        function BuildPressureResidual()
        end
        function BuildPressureMatrix()
        end
        function BuildTransportResidual()
        end
        function UpdatePressure()
        end
        function UpdateSaturation()
        end
    end
end