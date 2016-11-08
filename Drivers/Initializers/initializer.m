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
    end
    methods
        function InitializeProductionSystem(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            % Initialize Reservoir state           
            % 1. Initialize State object
            ProductionSystem.Reservoir.State.Initialize(DiscretizationModel.ReservoirGrid.N, FluidModel.NofPhases, FluidModel.NofComp);
            
            % 2. Assign Depth
            DiscretizationModel.ReservoirGrid.ComputeDepth(Formulation.GravityModel.alpha);
            
            % 3. Compute initial phase and component distribution
            ProductionSystem.Reservoir.State.T = ProductionSystem.Reservoir.Temp; % For now it's fine like this
            
            % 4. Define initial state
            Formulation.SinglePhase = obj.ComputeInitialState(ProductionSystem, FluidModel, Formulation.GravityModel, DiscretizationModel.ReservoirGrid);
            
            % Output initial status:      
            disp('Initial conditions:')
            disp(['reservoir pressure:' num2str(max(ProductionSystem.Reservoir.State.p/1e6)), ' MPa']);
            disp(['reservoir saturation:' num2str(max(ProductionSystem.Reservoir.State.S))]);
            disp(['reservoir z1: ', num2str(max(ProductionSystem.Reservoir.State.z(:,1)))])
            disp(['reservoir temperature: ', num2str(ProductionSystem.Reservoir.Temp)]);
            disp('---------------------------------------------------------');
            disp(char(5));
        end
        function SinglePhase = ComputeInitialState(obj, ProductionSystem, FluidModel, GravityModel, ReservoirGrid)
             % Define initial values
            P_init = ones(ReservoirGrid.N, 1)*1e5;
            z_init = ones(ReservoirGrid.N, 1)*0.037;
                        
            % 1. Assign initial valus
            ProductionSystem.Reservoir.State.p = P_init;
            ProductionSystem.Reservoir.State.z(:,1) = z_init;
            ProductionSystem.Reservoir.State.z(:,2) = 1 - z_init;
           
            % 2. Update Composition of the phases (Flash)
            SinglePhase = FluidModel.Flash(ProductionSystem.Reservoir.State);
            
            % 3 Compute Phase Density
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            
            % 4. Update S based on components mass balance
            FluidModel.ComputePhaseSaturation(ProductionSystem.Reservoir.State, SinglePhase);
            
            % 5. Total Density
            ProductionSystem.Reservoir.State.rhoT = FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State.S, ProductionSystem.Reservoir.State.rho);
            
            % 6. Compute initial Pc
            ProductionSystem.Reservoir.State.pc = FluidModel.ComputePc(ProductionSystem.Reservoir.State.S);
        end
    end
end
