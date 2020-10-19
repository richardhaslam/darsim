% Producer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef producer < handle
    properties
        Coordinate
        Cells
        p
        h
        qv % total volumetric rate
        QPhases
        QComponents
        Qh % enthalpy flux
        PI
        BHPDepth
    end
    methods
        function obj = producer(PI, coord)
            obj.Coordinate = coord;
            obj.PI = PI;
            obj.qv = 0;
        end
        function ResizeObjects(obj, n)
            obj.p =  ones(n,1) * obj.p;
            obj.qv = ones(n,1) * obj.qv/n;
        end
    end
end