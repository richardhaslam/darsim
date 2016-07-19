% Quadratic relative permeability model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
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
            % Non-Wetting phase relative permeability
            kr(:,2) = (1-S).^2;    
        end
        function dkr = ComputeDerivative(obj, Phases, s)
            S = (s-Phases(1).sr)/(1-Phases(1).sr-Phases(2).sr);
            dkr(:,1) = (1-Phases(1).sr-Phases(2).sr)^(-1)*2*S;
            dkr(:,2) = -(1-Fluid.sr(2)-Fluid.sr(1))^(-1)*2*(1-S);
        end
    end
end