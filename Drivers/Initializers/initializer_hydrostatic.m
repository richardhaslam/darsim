% Initializer base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 8 November 2016
%Last modified: 8 November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer_hydrostatic < initializer
    properties
    end
    methods
        function SinglePhase = ComputeInitialState(obj, ProductionSystem, FluidModel, GravityModel, ReservoirGrid)
             % Define initial values
            P_init = ones(ReservoirGrid.N, 1)*1e5;
            z_init = ones(ReservoirGrid.N, 1)*0.05;
            
            
            % Pressure
            P_init = P_init + Head;
            
            % 1. Assign initial valus
            ProductionSystem.Reservoir.State.p = P_init;
            ProductionSystem.Reservoir.State.z(:,1) = z_init;
            ProductionSystem.Reservoir.State.z(:,2) = 1 - z_init;
           
            SinglePhase = FluidModel.InitializeReservoir(ProductionSystem.Reservoir.State);
        end
    end
end

