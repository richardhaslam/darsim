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
        function k = Compute(obj, p, T, Components)
           k(:,1) = 1;
           k(:,2) = 0;
        end
        function dk = DKvalDp(obj, p)
            dk(:,1) =  1;
            dk(:,2) =  0;
        end
    end
end