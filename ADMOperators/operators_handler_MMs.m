% Operators handler for multilevel multiscale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 22 September 2017
%Last modified: 22 September 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef operators_handler_MMs < operators_handler
    properties
        R
        P
        C % Correction functions
    end
    methods
         function obj = operators_handler_MMs(cf)
            obj@operators_handler(cf);
        end
    end
end