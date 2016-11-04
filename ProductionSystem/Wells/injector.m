% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
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
        function obj = injector(PI, coord)
            obj.PI = PI;
            obj.Coord = coord;
        end
    end
end