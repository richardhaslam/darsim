% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 16 December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef injector < handle
    properties
        Coord
        Cells
        p
        T
        QPhases
        QComponents
        PI
        z
        S
        ni
        x
        x2
        rho
        Mob
    end
    methods
        function obj = injector(PI, coord, n_phases)
            obj.PI = PI;
            obj.Coord = coord;
            obj.rho = ones(1, n_phases);
        end
        function ResizeObjects(obj, n)
            obj.p =  ones(n,1) * obj.p;
            n_phases = length(obj.rho);
            obj.rho = ones(n, n_phases);
        end
    end
end