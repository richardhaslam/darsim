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
        Por0 % initial porosity
        Por % current porosity
        DPor
        Dp = 1e-3; % grain diameter
        K
        P0
        State
        State_old
        K_Cond_rock % rock conductivity
        Cpr % rock specific heat
        Cr % rock compressibility
        Rho % density of reservoir rock
        TotalPV
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
        function AddPermeabilityPorosity(obj, k, por0)
            obj.K = k;
            obj.Por0 = por0;
            obj.Por = por0;
            obj.TotalPV = obj.Length * obj.Width * obj.Thickness * obj.Por;
        end
    end
end