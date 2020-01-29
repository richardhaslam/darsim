% InitializerSinglePhase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer_Multiphase < initializer
    properties
    end
    methods
        function obj = initializer_Multiphase(names, values)
            obj@initializer(names, values);
        end
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started multi phase initialization');
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
                
            % DO WE NEED TO TAKE INTO ACCOUNT THE INITIAL BOUNDARY CONDITIONS, I.E. THE HEAT FLUX INTO THE RESERVOIR ??
            
            
            %% 3 Compute Phase Properties
            % Reservoir
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State); % call initial density
            switch FluidModel.name  
                case {'Geothermal_1T','Geothermal_2T'}
                FluidModel.ComputePhaseEnthalpies(ProductionSystem.Reservoir.State); % call enthalpy
                FluidModel.ComputePhaseViscosities(ProductionSystem.Reservoir.State); % call viscosity
                FluidModel.AddPhaseConductivities(ProductionSystem.Reservoir.State);
            end
            % Fractures
            if ProductionSystem.FracturesNetwork.Active
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    FluidModel.ComputePhaseDensities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    switch FluidModel.name
                        case {'Geothermal_1T','Geothermal_2T'}
                            FluidModel.ComputePhaseEnthalpies(ProductionSystem.FracturesNetwork.Fractures(f).State); % call enthalpy to fracture
                            FluidModel.ComputePhaseViscosities(ProductionSystem.FracturesNetwork.Fractures(f).State); % call viscosity to fracture
                            FluidModel.AddPhaseConductivities(ProductionSystem.FracturesNetwork.Fractures(f).State);
                    end
                end
            end
            
            % Output initial status:
            disp('Initial conditions:')
            
            % Initial status for reservoir:
            disp(['Pressure:' num2str(max(ProductionSystem.Reservoir.State.Properties('P_1').Value/1e5)), ' bar']);
            disp(['Saturation: ', num2str(1), ' (Single Phase)']);
            disp(['Temperature: ', num2str(ProductionSystem.Reservoir.Temp)]);
            disp('---------------------------------------------------------');
            
            disp(newline);
        end
    end
end