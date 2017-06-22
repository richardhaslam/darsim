% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 2 March 2017
%Last modified: 21 June 2017
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
            obj.BHPDepth = max(h);
            obj.p = obj.p - rho*GravityModel.g* (obj.BHPDepth - h);
        end
        function UpdateState(obj, State, K, FluidModel)
            p = State.Properties(['P_',num2str(FluidModel.NofPhases)]);
            for i = 1:FluidModel.NofPhases
                obj.QPhases(:,i) = obj.rho(:,i) .* obj.Mob(:,i) * obj.PI .* K(obj.Cells).* (obj.p - p.Value(obj.Cells));
            end
            obj.QComponents = zeros(length(obj.Cells), FluidModel.NofComp);
            switch(FluidModel.name)
                case('SinglePhase')
                case('Immiscible')
                otherwise
                    for j=1:FluidModel.NofComp
                        for phase=1:FluidModel.NofPhases
                            obj.QComponents(:, j) = obj.QComponents(:, j) + obj.x(:,(j-1)*2 + phase) .* obj.QPhases(:, phase);
                        end
                    end
            end
        end
        function [dQdp, dQdS] = dQPhasesdPdS(obj, State, K, NofPhases)
            p = State.Properties(['P_',num2str(NofPhases)]);
            dQdp = zeros(length(obj.Cells), NofPhases);
            dQdS = zeros(length(obj.Cells), NofPhases * (NofPhases - 1));
            for i = 1:NofPhases
                dQdp(:, i) = - obj.rho(:,i) .* obj.Mob(:,i) * obj.PI .* K(obj.Cells);
                dQdS(:, i) = 0;
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