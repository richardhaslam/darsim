% Producer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef producer_pressure < producer
    properties

    end
    methods
        function obj = producer_pressure(PI, coord, pressure)
            obj@producer(PI, coord)
            obj.p = pressure;
        end
        function UpdateState(obj, State, K, Mob, n_phases, n_components)
            for i = 1:n_phases
                obj.QPhases(:,i) = State.rho(obj.Cells, i) * Mob(obj.Cells, i) * obj.PI .* K(obj.Cells).* (obj.p - State.p(obj.Cells));
            end
            obj.QComponents(:, 1) = State.x1(obj.Cells, 1) * obj.QPhases(:,1) + State.x1(obj.Cells, 2) * obj.QPhases(:,2);
            obj.QComponents(:, 2) = (1 - State.x1(obj.Cells, 1)) * obj.QPhases(:,1) + (1 - State.x1(obj.Cells, 2)) * obj.QPhases(:,2);
        end
    end
end