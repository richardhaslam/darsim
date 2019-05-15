% Wilson K-values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 April 2017
%Last modified: 13 April 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Wilson_Kvalues_calculator < Kvalues_calculator
    properties
        T
    end
    methods
        function obj = Wilson_Kvalues_calculator(Temperature)
            obj.T = Temperature;
        end
        function K = Compute(obj, Status, Components, Phases)
            %% Compute K values
            P = Status.Properties('P_2').Value;
            N = length(P);
            n_c = length(Components);
            % Wilson correlation
            K = zeros(N, n_c);
            for i=1:n_c
                K(:, i) = Components(i).Pcrit./P * exp(5.37*(1+Components(i).w) .* (1 - Components(i).Tcrit./obj.T));
            end
        end
        function dK = DKvalDp(obj, Status, Components, Phases)
            %% Compute dKdp
            P = Status.Properties('P_2').Value;
            N = length(P);
            n_c = length(Components);
            % Derivative of Wilson correlation
            dK = zeros(N, n_c);
            for i=1:n_c
                dK(:, i) = - Components(i).Pcrit./P.^2 * exp(5.37*(1+Components(i).w) .* (1 - Components(i).Tcrit/obj.T));
            end
        end
    end
end
