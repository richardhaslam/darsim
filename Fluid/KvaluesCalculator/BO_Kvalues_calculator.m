% Black Oil K-values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 September 2016
%Last modified: 25 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_Kvalues_calculator < Kvalues_calculator
    properties
       Pref = 1e7;
    end
    methods
        function k = Compute(obj, Status, Components, Phases)
           [Rs, ~] = obj.ComputeRs(Status, Phases);
           k(:,1) = (Components(1).rho * Rs(:,2) + Components(2).rho * ones(length(Rs), 1)) ./ (Rs(:,2) * Components(1).rho);
           k(:,2) = 0;
        end
        function dk = DKvalDp(obj, Status, Components, Phases)
            [Rs, dRs] = obj.ComputeRs(Status, Phases);
            dk(:,1) =  - dRs(:,2) * Components(2).rho ./ (Components(1).rho * Rs(:,2).^2);
            dk(:,2) =  0;
        end
        function [Rs, dRs] = ComputeRs(obj, Status, Phases)
            n_phases = length(Phases);
            Pdim = Status.Properties(['P_',num2str(n_phases)]).Value/obj.Pref;
            Rs = zeros(length(Pdim), n_phases);
            dRs = zeros(length(Pdim), n_phases);
            for i=1:n_phases
                [Rs(:,i), dRs(:,i)] = Phases(i).ComputeRs(Pdim);
            end
            dRs = dRs / obj.Pref;
        end
    end
end