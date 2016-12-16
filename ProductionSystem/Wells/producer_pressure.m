% Producer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 16 December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef producer_pressure < producer
    properties

    end
    methods
        function obj = producer_pressure(PI, coord, pressure)
            obj@producer(PI, coord)
            obj.p = pressure;
        end
        function AdjustConstraint(obj, GravityModel, rhoT, h)
            rho = rhoT(obj.Cells);
            obj.p = obj.p - rho .*GravityModel.g .* h;
        end
        function UpdateState(obj, State, K, Mob, n_phases, n_components)
            for i = 1:n_phases
                obj.QPhases(:,i) = State.rho(obj.Cells, i) .* Mob(obj.Cells, i) * obj.PI .* K(obj.Cells).* (obj.p - State.p(obj.Cells));
            end
            for j=1:n_components
                obj.QComponents(:, j) = State.x(obj.Cells, (j-1)*2+1) .* obj.QPhases(:,1) + State.x(obj.Cells, (j-1)*2+2) .* obj.QPhases(:,2);
            end
        end
        function [A, rhs] = AddToPressureSystem(obj, Mob, K, A, rhs)
            a = obj.Cells;
            for ii=1:length(a)
                A(a(ii),a(ii)) = A(a(ii),a(ii)) + obj.PI * K(a(ii)) .* Mob(a(ii));
                rhs(a(ii)) = rhs(a(ii)) + obj.PI * K(a(ii)) .* Mob(a(ii)) .* obj.p;
            end
        end
        function q = TotalFlux(obj, q, p, K, Mob)
            q(obj.Cells) = q(obj.Cells) + obj.PI .* K(obj.Cells) .* Mob(obj.Cells,1) .* (obj.p - p(obj.Cells));
        end
    end
end