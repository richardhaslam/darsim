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
        K_Cond_eff % effective conductivity
        Cpr % rock specific heat
        rhoRock % density of reservoir rock
        Dp = 1e-3;% grain diameter
    end
    methods
        function obj = reservoir(length, width, thickness, temp)
            obj.Length = length;
            obj.Width = width;
            obj.Thickness = thickness;
            obj.Temp = temp;
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
        function AddConductivity(obj, k_cond_rock, k_cond_fluid)
            obj.K_Cond_rock = k_cond_rock;
            obj.K_Cond_eff  = k_cond_fluid * obj.Por + k_cond_rock * (1-obj.Por) * ones(size(obj.K));
        end
        function AddCoarsePermeability(obj, k_coarse)
            obj.K_coarse = k_coarse;
        end
    end
end