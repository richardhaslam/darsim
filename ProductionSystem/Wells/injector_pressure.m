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
            obj.p = obj.p - rho * GravityModel.g * (obj.BHPDepth - h);
        end
        function UpdateState(obj, State, K, FluidModel)
            for i = 1:FluidModel.NofPhases
                p = State.Properties(['P_',num2str(i)]).Value;
                obj.QPhases(:,i) = obj.rho(:,i) .* obj.Mob(:,i) .* obj.PI .* K(obj.Cells) .* (obj.p - p(obj.Cells));
            end
            obj.QComponents = zeros(length(obj.Cells), FluidModel.NofComponents);
            switch(FluidModel.name)
                case('SinglePhase')
                case('Immiscible')
                case{'Geothermal_SinglePhase','Geothermal_MultiPhase'}
                    for i = 1:FluidModel.NofPhases
                        obj.Qh(:,i) = obj.h(:,i) .* obj.QPhases(:,i);
                    end
                case('Compositional')
                    for j=1:FluidModel.NofComponents
                        for phase=1:FluidModel.NofPhases
                            obj.QComponents(:, j) = obj.QComponents(:, j) + obj.x(:,(j-1)*2 + phase) .* obj.QPhases(:, phase);
                        end
                    end
                otherwise
                    error('This fluid model is not implemented for injectors');
            end
        end
        %% Mass Flux Derivatives. IN THE INJECTOR, MOST PROPERTIES ARE CONSTANT (UNDER CONSTANT PRESSURE,TEMPERATURE or ENTHALPY)!!)
        function dQdp = ComputeWellMassFluxDerivativeWithRespectToPressure(obj, K, NofPhases)
            dQdp = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                dQdp(:, i) = obj.rho(:,i) .* obj.Mob(:,i) .* obj.PI .* K(obj.Cells) .* (-1);
            end
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
        end
        %% Heat Flux Derivatives. IN THE INJECTOR, MOST PROPERTIES ARE CONSTANT (UNDER CONSTANT PRESSURE,TEMPERATURE or ENTHALPY)!!)
        % Note for future developments: We are taking the derivatives in the wells, but using cell properties. We should use the well properties to take the derivatives
        function dQhdp = ComputeWellHeatFluxDerivativeWithRespectToPressure(obj, K, NofPhases)
            dQhdp = zeros(length(obj.Cells), NofPhases);
            for i = 1:NofPhases
                dQhdp(:, i) = obj.rho(:,i) .* obj.h(:,i) .* obj.PI .* K(obj.Cells) .* obj.Mob(:,i) .* (-1);  
            end
        end
        function dQhdS = ComputeWellHeatFluxDerivativeWithRespectToSaturation(obj, NofPhases)
            dQhdS = zeros(length(obj.Cells), NofPhases .* (NofPhases - 1));
            for i = 1:NofPhases
                dQhdS(:, i) = 0;
            end
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
        end
        %%
        function q = TotalFlux(obj, q, p, K)
            q(obj.Cells) = q(obj.Cells) + obj.PI .* K(obj.Cells) .* obj.Mob(:,1) .* (obj.p - p(obj.Cells));
        end
    end
end