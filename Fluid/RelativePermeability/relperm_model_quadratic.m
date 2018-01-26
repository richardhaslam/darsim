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
        n
    end
    methods
        function kr = ComputeRelPerm(obj, Phases, s)
            % Rescale saturations
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            S = max(S, 0);
            % Phase 1 relative permeability
            kr(:,1) = S.^obj.n(1);
            kr(s < Phases(1).sr, 1) = 0;
            kr(s < Phases(1).sr, 2) = 1;
            % Phase 2 relative permeability
            kr(:,2) = (1-S).^obj.n(2);
            kr(s > 1 - Phases(2).sr, 2) = 0;
            kr(s > 1 - Phases(2).sr, 1) = 1;
        end
        function dkr = ComputeDerivative(obj, Phases, s)
            S = zeros(length(s), 1);
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            
            dkr(:,1) = (1-Phases(1).sr-Phases(2).sr)^(-1)* obj.n * S.^(obj.n-1);
            dkr(s < Phases(1).sr, 1) = 0;
            dkr(s < Phases(1).sr, 2) = 0;
            dkr(:,2) = -(1-Phases(1).sr-Phases(2).sr)^(-1)*2*(1-S);
            dkr(s > 1 - Phases(2).sr, 2) = 0;
            dkr(s > 1 - Phases(2).sr, 1) = 0;
        end
        function ddkr = ComputeSecondDerivative(obj, Phases, s)
            ddkr(:,1) = ones(length(s), 1) * (1-Phases(1).sr-Phases(2).sr)^(-1)*2;
            ddkr(s < Phases(1).sr, 1) = 0;
            ddkr(s < Phases(1).sr, 2) = 0;
            ddkr(:,2) = ones(length(s), 1) * (1-Phases(1).sr-Phases(2).sr)^(-1)*2;
            ddkr(s > 1 - Phases(2).sr, 2) = 0;
            ddkr(s > 1 - Phases(2).sr, 1) = 0;
        end
    end
end