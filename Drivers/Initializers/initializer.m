% Initializer base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 8 November 2016
%Last modified: 8 November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer < handle
    properties
        VarNames
        VarValues
    end
    methods
        function obj = initializer(names, values)
            obj.VarNames = names;
            obj.VarValues = values;
        end
        function InitializeProductionSystem(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            % Initialize Reservoir state           
            % 1. Initialize State object (Has to be removed once the other one is ready)
            ProductionSystem.InitializeStates(DiscretizationModel, FluidModel.NofPhases, FluidModel.NofComp);
            
            % 1. Assign initial values 
            ProductionSystem.AssignInitialState(obj.VarNames, obj.VarValues);
            
            % 2. Compute initial phase and component distribution
            ProductionSystem.Reservoir.State.T = ProductionSystem.Reservoir.Temp; % For now it's fine like this
            
            % 3. Define initial state
            start = tic;
            obj.ComputeInitialState(ProductionSystem, FluidModel, Formulation, DiscretizationModel);
            initialization = toc(start);
            disp(['Initialization took ', num2str(initialization), ' s']);
            disp(char(5));
        end
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started simple initialization');
            
             % Define initial values
            P_init = ones(DiscretizationModel.ReservoirGrid.N, 1)*1e6;
            z_init = ones(DiscretizationModel.ReservoirGrid.N, 1)*0.037;
                        
            % 1. Assign initial values
            ProductionSystem.Reservoir.State.p = P_init;
            ProductionSystem.Reservoir.State.z(:,1) = z_init;
            ProductionSystem.Reservoir.State.z(:,2) = 1 - z_init;
           
            % 2. Update Composition of the phases (Flash)
            Formulation.SinglePhase = FluidModel.Flash(ProductionSystem.Reservoir.State);
            
            % 3 Compute Phase Density
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            
            % 4. Update S based on components mass balance
            FluidModel.ComputePhaseSaturation(ProductionSystem.Reservoir.State, Formulation.SinglePhase);
            
            % 5. Total Density
            ProductionSystem.Reservoir.State.rhoT = FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State.S, ProductionSystem.Reservoir.State.rho);
            
            % 6. Compute initial Pc
            ProductionSystem.Reservoir.State.pc = FluidModel.ComputePc(ProductionSystem.Reservoir.State.S);
            
            % Output initial status:      
            disp('Initial conditions:')
            disp(['reservoir pressure:' num2str(max(ProductionSystem.Reservoir.State.Properties('Pressure').Value/1e5)), ' bar']);
            disp(['reservoir saturation:' num2str(max(ProductionSystem.Reservoir.State.S))]);
            disp(['reservoir initial z: ', num2str(ProductionSystem.Reservoir.State.z(1,:))]);
            disp(['reservoir temperature: ', num2str(ProductionSystem.Reservoir.Temp)]);
            disp('---------------------------------------------------------');
            disp(char(5));
        end
    end
end

