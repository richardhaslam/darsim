% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 25 July 2016
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
        function UpdateState(obj, State, K, n_phases, n_components)
            for i = 1:n_phases
                obj.QPhases(:,i) = obj.rho(i) .* obj.Mob(:,i) * obj.PI .* K(obj.Cells).* (obj.p - State.p(obj.Cells));
            end
            for j=1:n_components
                obj.QComponents(:, j) = obj.x((j-1)*2 + 1) .* obj.QPhases(:,1) + obj.x((j-1)*2 + 2) .* obj.QPhases(:,2);
            end
        end
        function [A, rhs] = AddToPressureSystem(obj, K, A, rhs)
            a = obj.Cells;
            for ii=1:length(a)
                A(a(ii),a(ii)) = A(a(ii),a(ii)) + obj.PI * K(a(ii)) .* obj.Mob(ii, 1);
                rhs(a(ii)) = rhs(a(ii)) + obj.PI .* K(a(ii)) .* obj.Mob(ii,1) .* obj.p;
            end
        end
        function q = TotalFlux(obj, q, p, K)
            q(obj.Cells) = q(obj.Cells) + obj.PI .* K(obj.Cells) .* obj.Mob(:,1) .* (obj.p - p(obj.Cells));
        end
    end
end