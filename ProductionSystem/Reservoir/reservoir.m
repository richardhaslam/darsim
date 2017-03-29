% Reservoir 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 14 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reservoir < handle
    properties
        Length
        Width
        Thickness
        Temp
        K
        Por
        State
    end
    methods
        function obj = reservoir(length, width, thickness, temp)
            obj.Length = length;
            obj.Width = width;
            obj.Thickness = thickness;
            obj.Temp = temp;
            obj.State = status();
        end
        function AddPermeabilityPorosity(obj, k, por)
            obj.K = k;
            obj.Por = por;
        end
    end
end