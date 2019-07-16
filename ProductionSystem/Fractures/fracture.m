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
        k_cond
    end
    methods
        function obj = fracture()
            obj.State = status();
            obj.State_old = status();
        end
        function ComputePorosity(obj, P)
            obj.Por;
        end
        function ComputeDerPorosity(obj, P)
            obj.DPor = 0 * obj.Por;
        end
        function AddPermeabilityPorosity(obj, k, por)
            obj.K = k;
            obj.Por = por;
        end
    end
end