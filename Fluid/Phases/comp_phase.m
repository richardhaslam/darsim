% Phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 July 2016
%Last modified: 28 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef comp_phase < phase
    properties
        rho0 % Reference density
        cf % Compressibility
    end
    properties (Constant)
        Pref = 1e5; % Atmospheric pressure
    end
    methods
        function rho = ComputeDensity(obj, p, Components)
            rho = obj.rho0 .* exp(obj.cf.*(p - obj.Pref));
        end
        function drho = DrhoDp(obj, p, Components)
            drho = obj.cf .* obj.rho0 .*exp (obj.cf.*(p - obj.Pref));
        end
    end
end