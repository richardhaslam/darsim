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
        function UpdateState(obj, p, K, Mob)
            for i = 1:min(size(obj.QPhases))
                obj.QPhases(:,i) = Mob(obj.Cells, i) * obj.PI .* K(obj.Cells).* (obj.p - p(obj.Cells));
            end
        end
    end
end