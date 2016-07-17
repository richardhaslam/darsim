% Phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef phase < handle
    properties
        rho0 % Reference density
        mu % Reference Viscosity
        cf % Compressibility
        sr % Irriducible saturation
    end
    properties (Constant)
        Pref = 1e5; % Atmospheric pressure
    end
    methods
        function [rho, mu] = UpdatePhaseProperties(obj, Status)
            rho = ComputeDensity(obj, Status);
            mu = obj.mu;
        end
        function rho = ComputeDensity(obj, p)
            rho = obj.rho0 .* exp(obj.cf.*(p - obj.Pref));
        end
    end
end