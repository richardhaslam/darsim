% Linear relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_linear < relperm_model
    properties
    end
    methods
        function kr = ComputeRelPerm(obj, Phases, s)
            % Rescale saturations
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            kr(:,1) = S;
            kr(:,2) = 1-S;
        end
        function dkr = ComputeDerivative(obj, s)
            dkr(:,1) = (1-Phases(1).sr-Phases(2).sr)^(-1).*ones(length(s));
            dkr(:,2) = -(1-Fluid.sr(2)-Fluid.sr(1))^(-1).*ones(length(s));
        end
    end
end