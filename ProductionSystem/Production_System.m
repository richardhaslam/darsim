% Production System class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 3 March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Production_System < handle
    properties
        Reservoir
        FracturesNetwork = fracture_system();
        Wells
    end
    methods
        function AddReservoir(obj, Reservoir)
            obj.Reservoir = Reservoir;
        end
        function AddWells(obj, Wells)
            obj.Wells = Wells;
        end
        function AddFractures(obj, fractures)
            obj.FracturesNetwork = fractures;
        end
        function AssignInitialState(obj, VarNames, VarValues)
            obj.Reservoir.State.AssignInitialValues(VarNames, VarValues);
            if obj.FracturesNetwork.Active
                for f = 1:obj.FracturesNetwork.NumOfFrac
                    %% MODIFY INITIAL VALUES FOR FRAC
                    obj.FracturesNetwork.Fractures(f).State.AssignInitialValues(VarNames, VarValues);
                end
            end
        end
        function InitializeWells(obj, FluidModel, GravityModel, DiscretizationModel)
            %% Initialize Wells: 
            % 1. Create objects
            obj.Wells.InitializeFluxes(FluidModel.NofPhases, FluidModel.NofComp);
            % 2. Injectors pressures are adjusted to account for gravity
            % 2.a Estimate fluid properties inside injectors
            FluidModel.InitializeInjectors(obj.Wells.Inj);
            % Adjust pressures
            for i=1:obj.Wells.NofInj
                h = DiscretizationModel.ReservoirGrid.Depth(obj.Wells.Inj(i).Cells);
                obj.Wells.Inj(i).AdjustConstraint(GravityModel, h); 
            end
            % Create a gradient also inside the producers
            for i=1:length(obj.Wells.Prod)
                h = DiscretizationModel.ReservoirGrid.Depth(obj.Wells.Prod(i).Cells);
                obj.Wells.Prod(i).AdjustConstraint(GravityModel, obj.Reservoir.State.Properties('rhoT').Value, h);
            end
            
            % 3. Injection fluid properties are defined
            FluidModel.InitializeInjectors(obj.Wells.Inj);
            
            % 4. ComputeFluxes
            obj.Wells.UpdateState(obj.Reservoir, FluidModel);
        end
    end
end