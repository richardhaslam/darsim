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
        Depth
        Temp
        K
        Por
        State
    end
    methods
        function obj = reservoir(length, width, depth, temp)
            obj.Length = length;
            obj.Width = width;
            obj.Depth = depth;
            obj.Temp = temp;
            obj.State = status();
        end
        function AddPermeabilityPorosity(obj, k, por)
            obj.K = k;
            obj.Por = por;
        end
    end
end