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
            
            %% Initialize Wells: 
            % 1. Perforated cells
            DiscretizationModel.DefinePerforatedCells(obj.Wells);
            % 2. Create objects
            obj.Wells.InitializeFluxes(FluidModel.NofPhases, FluidModel.NofComp);
            % 3. Injection fluid properties are defined
            FluidModel.InitializeInjectors(obj.Wells.Inj);
            % 4. ComputeFluxes
            obj.Wells.UpdateState(obj.Reservoir, FluidModel);
        end
    end
end