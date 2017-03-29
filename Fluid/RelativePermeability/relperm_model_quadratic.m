% Quadratic relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef relperm_model_quadratic < relperm_model
    properties
    end
    methods
        function kr = ComputeRelPerm(obj, Phases, s)
            % Rescale saturations
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            % Wetting phase relative permeability
            kr(:,1) = S.^2;
            kr(s < Phases(1).sr, 1) = 0;
            kr(s < Phases(1).sr, 2) = 1;
            % Non-Wetting phase relative permeability
            kr(:,2) = (1-S).^2;
            kr(s > 1 - Phases(2).sr, 2) = 0;
            kr(s > 1 - Phases(2).sr, 1) = 1;
        end
        function dkr = ComputeDerivative(obj, Phases, s)
            S = zeros(length(s), 1);
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            dkr(:,1) = (1-Phases(1).sr-Phases(2).sr)^(-1)*2*S;
            dkr(s < Phases(1).sr, 1) = 0;
            dkr(s < Phases(1).sr, 2) = 0;
            dkr(:,2) = -(1-Phases(1).sr-Phases(2).sr)^(-1)*2*(1-S);
            dkr(s > 1 - Phases(2).sr, 2) = 0;
            dkr(s > 1 - Phases(2).sr, 1) = 0;
        end
    end
end