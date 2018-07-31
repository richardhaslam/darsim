% Thermal compressible phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo ...
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef therm_comp_phase < phase
    properties
        rho0 % Reference density
        cf % Compressibility
    end
    methods
        function rho = ComputeDensity(obj, p, T)
            % fill in whatever you like
        end
        function drho = DrhoDp(obj, p, T)
             % fill in whatever you like
        end
        function mu = ComputeViscosity()
        end
        function dmudT = DmuDT()
        end
        function dmudp = DmuDp()
        end
        function h = ComputeEnthalpy()
        end
        function dhdT = DhDT()
        end
        function dhdp = DhDp()
        end
    end
end