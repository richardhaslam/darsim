% Production System class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 16 December 2016
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
        function InitializeWells(obj, FluidModel, GravityModel, DiscretizationModel)
            %% Initialize Wells: 
            % 1. Create objects
            obj.Wells.InitializeFluxes(FluidModel.NofPhases, FluidModel.NofComp);
            % 2. Injectors pressures are adjusted to account for gravity
            % 2.a Estimate fluid properties inside injectors
            FluidModel.InitializeInjectors(obj.Wells.Inj);
            % Adjust pressures
            for i=1:length(obj.Wells.Inj)
                h = DiscretizationModel.ReservoirGrid.Depth(obj.Wells.Inj(i).Cells);
                obj.Wells.Inj(i).AdjustConstraint(GravityModel, h); 
            end
            % Create a gradient also inside the producers
            for i=1:length(obj.Wells.Prod)
                h = DiscretizationModel.ReservoirGrid.Depth(obj.Wells.Prod(i).Cells);
                obj.Wells.Prod(i).AdjustConstraint(GravityModel, obj.Reservoir.State.rhoT, h);
            end
            
            % 3. Injection fluid properties are defined
            FluidModel.InitializeInjectors(obj.Wells.Inj);
            
            % 4. ComputeFluxes
            obj.Wells.UpdateState(obj.Reservoir, FluidModel);
        end
    end
end