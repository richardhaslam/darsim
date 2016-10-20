% BO oil phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 29 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_oil_phase < phase
    properties
        Bo = 1;
    end
    methods
        function obj = BO_oil_phase()
            obj.mu = 1e-3;
            obj.sr = 0;
        end
        function rho = ComputeDensity(obj, Status, Components, Rs)
            rho = (Components(2).rho * ones(length(Rs),1) + Components(1).rho * Rs)/obj.Bo;
        end
        function drho = DrhoDp(obj, p, Components, dRs)
          drho = Components(1).rho * dRs / obj.Bo;  
        end
        function [Rs, dRs] = ComputeRs(obj, p)
            Rs = 0.2 * p + 0.2 * ones(length(p),1);
            dRs = 0.2;
        end
    end    
end