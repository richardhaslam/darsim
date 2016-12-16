% Producer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef producer < handle
    properties
        Coord
        Cells
        p
        QPhases
        QComponents
        PI
    end
    methods
        function obj = producer(PI, coord)
            obj.Coord = coord;
            obj.PI = PI;
        end
        function ResizeObjects(obj, n)
            obj.p =  ones(n,1) * obj.p;
        end
    end
end