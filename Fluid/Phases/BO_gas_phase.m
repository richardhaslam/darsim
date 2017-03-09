% BO gas phase class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 28 September 2016
%Last modified: 29 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef BO_gas_phase < phase
    properties
        Bg = 0.25;
    end
    methods
        function obj = BO_gas_phase()
            obj.mu = 1e-4;
            obj.sr = 0;
        end
        function rho = ComputeDensity(obj, p, Components, rs)
            N = length(p);
            rho = Components(1).rho/obj.Bg * ones(N,1);
        end
        function drho = DrhoDp(obj, Status, Components, drs)
            drho = drs;
        end
        function [Rs, dRs] = ComputeRs(obj, p)
            Rs = 0;
            dRs = 0;
        end
    end    
end