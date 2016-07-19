% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef injector_pressure < injector
    properties
    end
    methods
        function obj = injector_pressure(PI, coord, pressure, temperature)
            obj@injector(PI, coord)
            obj.p = pressure;
            obj.T = temperature;
        end
        function UpdateState(obj, p, K)
            for i = 1:min(size(obj.QPhases))
                obj.QPhases(:,i) = obj.Mob(:,i) * obj.PI .* K(obj.Cells).* (obj.p - p(obj.Cells));
            end
        end
    end
end