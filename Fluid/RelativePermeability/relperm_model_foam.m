% Foam relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Esat Unal
%TU Delft
%Created: 25 September 2017
%Last modified: 25 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_foam < relperm_model
    properties
        fmob
        fmdry 
    end
    methods
        function kr = ComputeRelPerm(obj, Phases, s)
         
        end
        function dkr = ComputeDerivative(obj, Phases, s)
            
        end
        function ddkr = ComputeSecondDerivative(obj, Phases, s)
           
        end
    end
end