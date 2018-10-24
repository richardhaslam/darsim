% Producer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
            obj.BHPDepth = max(h);
            obj.p = obj.p - rho .*GravityModel.g .*(obj.BHPDepth - h);
        end
        function UpdateState(obj, State, K, Mob, FluidModel)
            p = State.Properties(['P_',num2str(FluidModel.NofPhases)]);
            for i = 1:FluidModel.NofPhases
                rho = State.Properties(['rho_', num2str(i)]);
                obj.QPhases(:,i) = rho.Value(obj.Cells) .* Mob(obj.Cells, i) * obj.PI .* K(obj.Cells).* (obj.p - p.Value(obj.Cells));
            end
            obj.QComponents = zeros(length(obj.Cells), FluidModel.NofComp);
            switch(FluidModel.name)
                case('SinglePhase')
                case('Immiscible')
                case('Geothermal')
                    % update also Qh flux
                otherwise
                    for j=1:FluidModel.NofComp
                        for phase=1:FluidModel.NofPhases
                            x = State.Properties(['x_', num2str(j), 'ph', num2str(phase)]);
                            obj.QComponents(:, j) = obj.QComponents(:, j) + x.Value(obj.Cells) .* obj.QPhases(:,phase);
                        end
                    end
            end
        end
        function [dQdp, dQdS] = dQPhasesdPdS(obj, State, K, Mob, dMob, drho, NofPhases)
            p = State.Properties(['P_',num2str(NofPhases)]);
            dQdp = zeros(length(obj.Cells), NofPhases);
            dQdS = zeros(length(obj.Cells), NofPhases * (NofPhases - 1));
            for i = 1:NofPhases
                rho = State.Properties(['rho_', num2str(i)]);
                dQdp(:, i) = - rho.Value(obj.Cells) .* Mob(obj.Cells,i) * obj.PI .* K(obj.Cells) + drho(obj.Cells, i) .* Mob(obj.Cells, i) * obj.PI .* K(obj.Cells).* (obj.p - p.Value(obj.Cells));
                dQdS(:, i) = rho.Value(obj.Cells) .* dMob(obj.Cells, i) * obj.PI .* K(obj.Cells).* (obj.p - p.Value(obj.Cells));
            end
        end
        function [dQdp, dQdT] = dQdPdT(obj, K, NofPhases) % need perforated cell properties
        end
        function [dQhdp, dQhdT] = dQhdPdT(obj, K, NofPhases) % need perforated cell properties
        end
%         function [A, rhs] = AddToPressureSystem(obj, Mob, K, A, rhs)
%             a = obj.Cells;
%             for ii=1:length(a)
%                 A(a(ii),a(ii)) = A(a(ii),a(ii)) + obj.PI * K(a(ii)) .* Mob(a(ii));
%                 rhs(a(ii)) = rhs(a(ii)) + obj.PI * K(a(ii)) .* Mob(a(ii)) .* obj.p;
%             end
%         end
        function q = TotalFlux(obj, q, p, K, Mob)
            q(obj.Cells) = q(obj.Cells) + obj.PI .* K(obj.Cells) .* Mob(obj.Cells,1) .* (obj.p - p(obj.Cells));
        end
    end
end