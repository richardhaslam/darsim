% Black Oil K-values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 September 2016
%Last modified: 28 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_Kvalues_calculator < Kvalues_calculator
    properties
       
    end
    methods
        function k = Compute(obj, Components, Rs)
           k(:,1) = (Components(1).rho * Rs(:,2) + Components(2).rho * ones(length(Rs), 1)) ./ (Rs(:,2) * Components(1).rho);
           k(:,2) = 0;
        end
        function dk = DKvalDp(obj, Components, Rs, dRs)
            dk(:,1) =  - dRs(:,2) * Components(2).rho ./ (Components(1).rho * Rs(:,2).^2);
            dk(:,2) =  0;
        end
    end
end