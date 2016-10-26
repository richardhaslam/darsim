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
        function k = Compute(obj, Status, Components, Phases)
           k(:,1) = 1.5 * ones(length(Status.p), 1);
           k(:,2) = 0.5 * ones(length(Status.p), 1);
        end
        function dk = DKvalDp(obj, Status, Components, Phases)
            dk(:,1) =  zeros(length(Status.p), 1);
            dk(:,2) =  zeros(length(Status.p), 1);
        end
    end
end