% BO gas phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 September 2016
%Last modified: 16 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_gas_phase < phase
    properties
        Pb = 2.5e7;
        Bg = 0.0025;
        Pst = 1e5;
        alpha
    end
    methods
        function obj = BO_gas_phase()
            obj.mu = 1e-4;
            obj.sr = 0; % it's assigned later on
            obj.alpha = log(obj.Bg)/(obj.Pb/obj.Pst - 1)*1/obj.Pst; 
        end
        function rho = ComputeDensity(obj, p, Components, rs)
            %N = length(p);
            %Bg_r = 0.25 * ones(length(Bg_r),1);
            %Bg_r = exp(obj.alpha .* (p - obj.Pst));
            Bg_r = 0.5*(p./obj.Pst).^-1;  
            rho = Components(1).rho./Bg_r;
        end
        function drho = ComputeDrhoDp(obj, p, Components, rs, drs)
            %Bg_r = exp(obj.alpha .* (p - obj.Pst));
            %Num = -obj.alpha * Bg_r;
            %drho = zeros(length(Bg_r), 1);
            Bg_r = 0.5*(p./obj.Pst).^-1;
            Num = 0.5*obj.Pst * p.^-2;
            Den = (Bg_r).^2;
            drho = Num ./ Den;
        end
        function [Rs, dRs] = ComputeRs(obj, p)
            Rs = 0;
            dRs = 0;
        end
        function [Rs, dRs] = RsOfUnderSaturatedPhase(obj, z, Components, Rs, dRs, SinglePhase)
            
        end
    end    
end