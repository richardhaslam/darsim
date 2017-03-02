% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 2 March 2017
%Last modified: 2 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef injector_pressure < injector
    properties
    end
    methods
        function obj = injector_pressure(PI, coord, pressure, temperature, n_phases)
            obj@injector(PI, coord, n_phases)
            obj.p = pressure;
            obj.T = temperature;
        end
        function AdjustConstraint(obj, GravityModel, h)
            rho = max(max(obj.rho));
            obj.p = obj.p - rho*GravityModel.g*h;
        end
        function UpdateState(obj, State, K, n_phases, n_components)
            for i = 1:n_phases
                obj.QPhases(:,i) = obj.rho(:,i) .* obj.Mob(:,i) * obj.PI .* K(obj.Cells).* (obj.p - State.p(obj.Cells) + (i - 3) * (1-i) * State.pc(obj.Cells));
            end
            obj.QComponents = obj.QComponents*0;
            for j=1:n_components
                for phase =1:n_phases
                    obj.QComponents(:, j) = obj.QComponents(:, j) + obj.x((j-1)*2 + phase) .* obj.QPhases(:,phase); 
                end
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