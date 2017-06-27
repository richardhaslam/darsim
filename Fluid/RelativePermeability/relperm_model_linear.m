% Linear relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_linear < relperm_model
    properties
    end
    methods
        function kr = ComputeRelPerm(obj, Phases, s)
            % Rescale saturations
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            kr(:,1) = S;
            kr(s < Phases(1).sr, 1) = 0;
            kr(:,2) = 1-S;
            kr(s > 1 - Phases(2).sr, 2) = 0; 
        end
        function dkr = ComputeDerivative(obj, Phases, s)
            dkr(:,1) = (1-Phases(1).sr - Phases(2).sr)^(-1) .* ones(length(s),1);
            dkr(s < Phases(1).sr, 1) = 0;
            dkr(:,2) = -(1-Phases(1).sr - Phases(2).sr)^(-1) .* ones(length(s),1);
            dkr(s > 1 - Phases(2).sr, 2) = 0;
        end
        function dkr = ComputeSecondDerivative(obj, Phases, s)
            dkr(:,1) = zeros(length(s), 1);
            dkr(s < Phases(1).sr, 1) = 0;
            dkr(:,2) = zeros(length(s), 1);
            dkr(s > 1 - Phases(2).sr, 2) = 0;
        end
    end
end