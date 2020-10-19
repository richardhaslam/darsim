% Fracture system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fracture < handle
    properties
        Length
        Width
        Thickness
        Temp
        Por
        DPor
        Dp = 1e-3; % grain diameter
        K
        State
        State_old
        K_Cond_rock % rock conductivity
        K_Cond_eff % effective conductivity
        Cpr % rock specific heat
        Rho % density of reservoir rock
    end
    methods
        function obj = fracture()
            obj.State = status();
            obj.State_old = status();
        end
        function ComputePorosity(obj, P)
            obj.Por = obj.Por0.*exp(obj.Cr.*(P-obj.P0));
        end
        function ComputeDerPorosity(obj, P)
            obj.DPor = obj.Cr.*obj.Por0.*exp(obj.Cr.*(P-obj.P0));
        end
        function AddPermeabilityPorosity(obj, k, por)
            obj.K = k;
            obj.Por = por;
        end
    end
end