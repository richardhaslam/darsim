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
        qv % total volumetric rate
        QPhases
        QComponents
        PI
        BHPDepth
    end
    methods
        function obj = producer(PI, coord)
            obj.Coord = coord;
            obj.PI = PI;
            obj.qv = 0;
        end
        function ResizeObjects(obj, n)
            obj.p =  ones(n,1) * obj.p;
            obj.qv = ones(n,1) * obj.qv/n;
        end
    end
end