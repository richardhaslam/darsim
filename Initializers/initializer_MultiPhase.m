% InitializerSinglePhase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer_MultiPhase < initializer
    properties
    end
    methods
        function obj = initializer_MultiPhase(names, values)
            obj@initializer(names, values);
        end
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started MultiPhase initialization');
            
            if isnan(ProductionSystem.Reservoir.State.Properties('hTfluid').Value)
                P = ProductionSystem.Reservoir.State.Properties('P_2').Value;
                T = ProductionSystem.Reservoir.State.Properties('T').Value;
                H = ProductionSystem.Reservoir.State.Properties('hTfluid');
                H.Value = FluidModel.Phases(1).ComputeWaterEnthalpy(P, T);
            end
            
            Formulation.ComputeProperties(ProductionSystem, FluidModel);
            
            % Output initial status:
            disp('Initial conditions:')
            
            % Initial status for reservoir:
            disp(['Pressure:' num2str(max(ProductionSystem.Reservoir.State.Properties('P_1').Value/1e5)), ' bar']);
            disp(['Enthalpy:']);
            disp(['Saturation: ', num2str(1), ' (Single Phase)']); % saturation of water and steam
            disp(['Temperature: ', num2str(ProductionSystem.Reservoir.Temp)]);
            disp('---------------------------------------------------------');
            
            disp(newline);
        end
    end
end