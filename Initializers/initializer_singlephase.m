% InitializerSinglePhase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer_singlephase < initializer
    properties
    end
    methods
        function obj = initializer_singlephase(names, values)
            obj@initializer(names, values);
        end
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started single phase initialization');
            %% 1. Change initial 
            
            %% 2. Update Composition of the phases (Flash)
            % Reservoir
            Formulation.SinglePhase(1) = FluidModel.Flash(ProductionSystem.Reservoir.State);
            % Fractures
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    Formulation.SinglePhase(f+1) = FluidModel.Flash(ProductionSystem.FracturesNetwork.Fractures(f).State);
                end
            end
                
            %% 3 Compute Phase Density
            % Reservoir
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            % Fractures
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                end
            end
            
            % Output initial status:
            disp('Initial conditions:')
            
            % Initial status for reservoir:
            disp(['reservoir pressure:' num2str(max(ProductionSystem.Reservoir.State.Properties('P_1').Value/1e5)), ' bar']);
            disp(['Single phase reservoir']);
            disp(['reservoir temperature: ', num2str(ProductionSystem.Reservoir.Temp)]);
            disp('---------------------------------------------------------');
            
            % Initial status for fractures:
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    disp(['Fracture ', num2str(f), ' pressure:', num2str(max(ProductionSystem.FracturesNetwork.Fractures(f).State.Properties('P_1').Value/1e5)), ' bar']);
                    disp(['Single phase fracture']);
                    disp(['Fracture ', num2str(f), ' temperature: ', num2str(ProductionSystem.FracturesNetwork.Fractures(f).Temp)]);
                    disp('---------------------------------------------------------');
                end
            end
            
            disp(char(5));
        end
    end
end