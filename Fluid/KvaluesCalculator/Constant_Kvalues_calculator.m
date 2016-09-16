% Standing K-values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 16 September 2016
%Last modified: 16 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Constant_Kvalues_calculator < Kvalues_calculator
    properties
       
    end
    methods
        function k = Compute(obj, p, T, Components)
           k(:,1) = 1.5 * ones(length(p), 1);
           k(:,2) = .5 * ones(length(p), 1);
        end
        function dk = DKvalDp(obj, p)
            dk(:,1) =  zeros(length(p), 1);
            dk(:,2) =  zeros(length(p), 1);
        end
    end
end