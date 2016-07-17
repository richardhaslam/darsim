% Production System class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 17 July 2016
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
            obj.Reservoir.State = FluidModel.InitializeReservoir(obj.Reservoir.State);
            
            %% Initialize Wells: Injection fluid properties have to be defined
            obj.Wells.Inj = FluidModel.InitializeInjectors(obj.Wells.Inj);
        end
    end
end