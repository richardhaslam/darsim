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
        function InitializeWells(obj, FluidModel)
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