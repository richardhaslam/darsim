% BO oil phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 15 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_oil_phase < phase
    properties
        Pb = 2.5e7;
        Pst = 1e5;
        m
        q
    end
    methods
        function obj = BO_oil_phase()
            obj.mu = 1e-3;
            obj.sr = 0;
            obj.m = 0.2/(obj.Pb/obj.Pst - 1) * 1/obj.Pst;
            obj.q = 1 - obj.m;
        end
        function rho = ComputeDensity(obj, p, Components, Rs)
            Bo_r = 1;
            %Bo_r = obj.m .* p + obj.q;
            rho = (Components(2).rho + Components(1).rho .* Rs)./Bo_r;
        end
        function drho = ComputeDrhoDp(obj, p, Components, Rs, dRs)
            %Bo_r = obj.m .* p + obj.q;
            %dBo_r = obj.m;
            Bo_r = 1;
            dBo_r = 0;
            Num = Bo_r .* (Components(1).rho * dRs) - dBo_r .* (Components(2).rho + Components(1).rho .* Rs);
            Den = Bo_r.^2;
            drho = Num ./ Den;
        end
        function [Rs, dRs] = ComputeRs(obj, p)
            Rs = 100 * p/obj.Pb + 0.2;
            dRs = 100 / obj.Pb;
        end
        function drho = ComputeDrhoDz(obj, p, z, Components, SinglePhase)
            %Bo_r = obj.m .* p + obj.q;
            Bo_r = ones(length(p), 1);
            drho = zeros(length(p), 1);
            dRs = zeros(length(p), 1);
            dRs(SinglePhase == 2) = Components(2).rho* Components(1).rho ./ ((1-z(SinglePhase == 2))*Components(1).rho).^2;
            drho(SinglePhase == 2) = Components(1).rho .* dRs((SinglePhase == 2)) ./ Bo_r (SinglePhase == 2);
        end
        function [Rs, dRs] = RsOfUnderSaturatedPhase(obj, z, Components, Rs, dRs, SinglePhase)
            Rs(SinglePhase == 2) = Components(2).rho .* z(SinglePhase == 2) ./ ((1 - z(SinglePhase == 2)) .* Components(1).rho);
            dRs(SinglePhase == 2) = 0;
        end
    end    
end