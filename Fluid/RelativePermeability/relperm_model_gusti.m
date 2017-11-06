% Gusti relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: All DARSim2 members
%TU Delft
%Created: 6 November 2017
%Last modified: 6 November 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_gusti < relperm_model
    properties
        a
    end
    methods
        function obj = relperm_model_gusti(input)
            obj.a = input;
        end
        function kr = ComputeRelPerm(obj, Phases, s)
            % Rescale saturations
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            S = max(S, 0);
            kr(:,1) = obj.a * S;
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