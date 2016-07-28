% Production System class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Production_System < handle
    properties
        Reservoir
        Wells
    end
    methods
        function AddReservoir(obj, Reservoir)
            obj.Reservoir = Reservoir;
        end
        function AddWells(obj, Wells)
            obj.Wells = Wells;
        end
        function Initialize(obj, DiscretizationModel, FluidModel)
            %% Initialize Reservoir state           
            % 1. Initialize State object
            obj.Reservoir.State.Initialize(DiscretizationModel.ReservoirGrid.N);
            % 2. Compute initial phase and component distribution
            obj.Reservoir.State.T = obj.Reservoir.Temp; % For now it's fine like this
            FluidModel.InitializeReservoir(obj.Reservoir.State);
            
            % Output initial status:      
            disp('Initial conditions:')
            disp(['reservoir pressure:' num2str(max(obj.Reservoir.State.p/1e6)), ' MPa']);
            disp(['reservoir saturation:' num2str(max(obj.Reservoir.State.S))]);
            disp(['reservoir z1: ', num2str(max(obj.Reservoir.State.z(:,1)))])
            disp(['reservoir temperature: ', num2str(obj.Reservoir.Temp)]);
            disp('---------------------------------------------------------');
            disp(char(5));
            
            %% Initialize Wells: 
            % 1. Create objects
            obj.Wells.InitializeFluxes(FluidModel.NofPhases, FluidModel.NofComp);
            % 2. Injection fluid properties are defined
            FluidModel.InitializeInjectors(obj.Wells.Inj);
            % 3. ComputeFluxes
            obj.Wells.UpdateState(obj.Reservoir, FluidModel);
        end
    end
end