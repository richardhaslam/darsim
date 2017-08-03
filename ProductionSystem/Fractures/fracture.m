% Fracture system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 July 2016
%Last modified: 03 August 2017
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
        PointA
        PointB
        PointC
        PointD
        GridCoords
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