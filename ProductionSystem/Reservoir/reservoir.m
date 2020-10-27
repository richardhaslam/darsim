% Reservoir 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reservoir < handle
    properties
        Length
        Width
        Thickness
        Temp
        K
        K_coarse
        Por0 % initial porosity
        Por % current porosity
        DPor
        TotalPV
        State
        level
        State_old
        MaxLevel
        P0
        Cr % rock compressibility
        K_Cond_rock % rock conductivity
        Cpr % rock specific heat
        Rho % density of reservoir rock
        Dp = 1e-3;% grain diameter
    end
    methods
        function obj = reservoir(length, width, thickness, temperature)
            obj.Length = length;
            obj.Width = width;
            obj.Thickness = thickness;
            obj.Temp = temperature;
            obj.State = status();
            obj.State_old = status();
        end
        function ComputePorosity(obj, P)
            obj.Por = obj.Por0.*exp(obj.Cr.*(P-obj.P0));
        end
        function ComputeDerPorosity(obj, P)
            obj.DPor = obj.Cr.*obj.Por0.*exp(obj.Cr.*(P-obj.P0));
        end
        function AddPermeabilityPorosity(obj, k, por0)
            obj.K = k;
            obj.Por0 = por0;
            obj.Por = por0;
            obj.TotalPV = obj.Length * obj.Width * obj.Thickness * obj.Por;
        end
        function AddCoarsePermeability(obj, k_coarse)
            obj.K_coarse = k_coarse;
        end
    end
end