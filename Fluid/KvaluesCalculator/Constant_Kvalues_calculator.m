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
           k = zeros(length(Status.p), length(Components));
           for c=1:length(Components)
               k(:, c) = Components(c).kval0 * ones(length(Status.p), 1);
           end
        end
        function dk = DKvalDp(obj, Status, Components, Phases)
            dk = zeros(length(Status.p), length(Components));
        end
    end
end