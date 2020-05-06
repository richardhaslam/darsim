% Injector 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
            for i = 1:FluidModel.NofPhases
                p = State.Properties(['P_',num2str(i)]).Value;
                obj.QPhases(:,i) = obj.rho(:,i) .* obj.Mob(:,i) .* obj.PI .* K(obj.Cells).* (obj.p - p(obj.Cells));
            end
            obj.QComponents = zeros(length(obj.Cells), FluidModel.NofComp);
            switch(FluidModel.name)
                case('SinglePhase')
                case('Immiscible')
                case{'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                    for i = 1:FluidModel.NofPhases
                        obj.Qh(:,i) = obj.h(:,i) .* obj.QPhases(:,i);
                    end
                otherwise
                    for j=1:FluidModel.NofComp
                        for phase=1:FluidModel.NofPhases
                            obj.QComponents(:, j) = obj.QComponents(:, j) + obj.x(:,(j-1)*2 + phase) .* obj.QPhases(:, phase);
                        end
                    end
            end
        end
        %% IN THE INJECTOR, MOST PROPERTIES ARE CONSTANT (UNDER CONSTANT PRESSURE,ENTHALPY) !!
        
        % Mass Flux Derivatives
        function dQdp = ComputeWellMassFluxDerivativeWithRespectToPressure(obj, K, NofPhases)
            dQdp = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                dQdp(:, i) = obj.rho(:,i) .* obj.Mob(:,i) .* obj.PI .* K(obj.Cells) .* (-1);
            end
%             dQdp(:, i) = drhodp(obj.Cells,i) .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells)) + ...
%                          rho(obj.Cells)      .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (-1                  ); % + ...
% %                             rho(obj.Cells)      .* obj.PI .* K(obj.Cells) .* dMobdp(obj.Cells,i) .* (obj.p - p(obj.Cells));
        end
        function dQdS = ComputeWellMassFluxDerivativeWithRespectToSaturation(obj, NofPhases)
            dQdS = zeros(length(obj.Cells), NofPhases .* (NofPhases - 1));
            for i = 1:NofPhases
                dQdS(:, i) = 0;
            end
        end
        function dQdT = ComputeWellMassFluxDerivativeWithRespectToTemperature(obj, NofPhases)
            dQdT = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                dQdT(:, i) = 0;   
            end
        end
        function dQdh = ComputeWellMassFluxDerivativeWithRespectToEnthalpy(obj, NofPhases)
            dQdh = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                dQdh(:, i) = 0;   
            end
%             dQdh(:, i) = 0; %drhodh(obj.Cells,i) .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells)); % + ...
% %                              rho(obj.Cells)      .* obj.PI .* K(obj.Cells) .* dMobdh(obj.Cells,i) .* (obj.p - p(obj.Cells));
        end
        
        % Heat Flux Derivatives
        % We are taking the derivatives in the wells, but using cell properties. We should use the well properties to take the derivatives
        function dQhdp = ComputeWellHeatFluxDerivativeWithRespectToPressure(obj, K, NofPhases)
            dQhdp = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                dQhdp(:, i) = obj.rho(:,i) .* obj.h(:,i) .* obj.PI .* K(obj.Cells) .* obj.Mob(:,i) .* (-1) ;  
            end
%             dQhdp(:, i) = drho_times_hdp(obj.Cells,i)    .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells)) + ...
%                           rho(obj.Cells) .* h(obj.Cells) .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (-1                  ) + ...
%                               rho(obj.Cells) .* h(obj.Cells) .* obj.PI .* K(obj.Cells) .* dMobdp(obj.Cells,i) .* (obj.p - p(obj.Cells));
        end
        function dQhdT = ComputeWellHeatFluxDerivativeWithRespectToTemperature(obj, NofPhases)
            dQhdT = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                dQhdT(:, i) = 0;
            end
        end
        function dQhdh = ComputeWellHeatFluxDerivativeWithRespectToEnthalpy(obj, NofPhases)
            dQhdh = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                dQhdh(:, i) = 0; 
            end
%             dQhdh(:, i) = drhodh(obj.Cells,i) .* h(obj.Cells) .* obj.PI .* K(obj.Cells) .* Mob(obj.Cells,i)    .* (obj.p - p(obj.Cells)) + ...
%                               rho(obj.Cells)    .* h(obj.Cells) .* obj.PI .* K(obj.Cells) .* dMobdh(obj.Cells,i) .* (obj.p - p(obj.Cells));
        end
        
        function q = TotalFlux(obj, q, p, K)
            q(obj.Cells) = q(obj.Cells) + obj.PI .* K(obj.Cells) .* obj.Mob(:,1) .* (obj.p - p(obj.Cells));
        end
    end
end