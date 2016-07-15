% Reservoir 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 8 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef reservoir < handle
    properties
        Length
        Width
        Depth
        Temp
        K
        por
    end
    methods
        function obj = reservoir(length, width, depth, temp)
            obj.Length = length;
            obj.Width = width;
            obj.Depth = depth;
            obj.Temp = temp;
        end
    end
end