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
        K
        State
        State_old
    end
    methods
        function obj = fracture()
            obj.State = status();
            obj.State_old = status();
        end
        function AddPermeabilityPorosity(obj, k, por)
            obj.K = k;
            obj.Por = por;
        end
    end
end