% Producer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef producer_rate < producer
    properties
    end
    methods
          function obj = producer_rate(PI, coord, rate)
            obj@producer(PI, coord)
            obj.qv = rate;
        end
        function AdjustConstraint(obj, GravityModel, rhoT, h)
           
        end
        function UpdateState(obj, State, K, Mob, FluidModel)        
            p = State.Properties(['P_',num2str(FluidModel.NofPhases)]); 
            obj.p = obj.qv ./ (obj.PI .* sum(Mob(obj.Cells), 2) .* K(obj.Cells)) + p.Value(obj.Cells);
            f = Mob(obj.Cells, :) ./ sum(Mob(obj.Cells,2), 2);
            for i = 1:FluidModel.NofPhases
                rho = State.Properties(['rho_', num2str(i)]);
                obj.QPhases(:,i) = rho(obj.Cells, i) .* f(:, i) * obj.qv; 
            end
            
            obj.QComponents = zeros(length(obj.Cells), FluidModel.NofComp);
%             switch(FluidModel.name)
%                 case('SinglePhase')
%                 case('Immiscible')
%                 otherwise
%                     for j=1:FluidModel.NofComp
%                         for phase=1:FluidModel.NofPhases
%                             x = State.Properties(['x_', num2str(j), 'ph', num2str(phase)]);
%                             obj.QComponents(:, j) = obj.QComponents(:, j) + x.Value(obj.Cells) .* obj.QPhases(:,phase);
%                         end
%                     end
%             end
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
        function [A, rhs] = AddToPressureSystem(obj, Mob, K, A, rhs)
            a = obj.Cells;
            for ii=1:length(a)
                A(a(ii),a(ii)) = A(a(ii),a(ii)) + obj.PI * K(a(ii)) .* Mob(a(ii));
                rhs(a(ii)) = rhs(a(ii)) + obj.PI * K(a(ii)) .* Mob(a(ii)) .* obj.p;
            end
        end
        function Q = TotalFlux(obj, Q, p, K, Mob)
            Q(obj.Cells) = Q(obj.Cells) + obj.PI .* K(obj.Cells) .* Mob(obj.Cells,1) .* (obj.p - p(obj.Cells));
        end  
    end
end