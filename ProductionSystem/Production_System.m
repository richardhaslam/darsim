% Production System class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
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
            obj.Reservoir.State_old.AssignInitialValues(VarNames, VarValues);
            if obj.FracturesNetwork.Active
                for f = 1:obj.FracturesNetwork.NumOfFrac
                    %% MODIFY INITIAL VALUES FOR FRAC
                    
                    obj.FracturesNetwork.Fractures(f).State.AssignInitialValues(VarNames, VarValues(1,:));
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
                depth = DiscretizationModel.ReservoirGrid.Depth(obj.Wells.Inj(i).Cells);
                obj.Wells.Inj(i).AdjustConstraint(GravityModel, depth); 
            end
            % Create a gradient also inside the producers
            for i=1:length(obj.Wells.Prod)
                depth = DiscretizationModel.ReservoirGrid.Depth(obj.Wells.Prod(i).Cells);
                obj.Wells.Prod(i).AdjustConstraint(GravityModel, obj.Reservoir.State.Properties('rhoT').Value, depth);
            end
            
            % 3. Injection fluid properties are defined
            FluidModel.InitializeInjectors(obj.Wells.Inj);
            
            % 4. ComputeFluxes
            obj.Wells.UpdateState(obj.Reservoir, FluidModel);
        end
        function SavePreviousState(obj)
            %% Reservoir
            obj.Reservoir.State_old.CopyProperties(obj.Reservoir.State);
            %% Fractures
            if obj.FracturesNetwork.Active
                for f = 1:obj.FracturesNetwork.NumOfFrac
                    %% MODIFY INITIAL VALUES FOR FRAC
                    obj.FracturesNetwork.Fractures(f).State_old.CopyProperties(obj.FracturesNetwork.Fractures(f).State);
                end
            end
        end
        function x_Global = CreateGlobalVariables(obj, FineGrid, n_phases, key)
            %
            x_Global   = zeros(sum([FineGrid(1:end).N]), n_phases);
            %% Reservoir
            for i=1:n_phases
                x = obj.Reservoir.State.Properties([key, num2str(i)]).Value;
                x_Global(1:FineGrid(1).N, i) = x;
            end
            %% Fractures
            if obj.FracturesNetwork.Active
                End = FineGrid(1).N;
                for f=1:obj.FracturesNetwork.NumOfFrac
                    Start = End+1;
                    End = Start + FineGrid(1+f).N - 1;
                    for i=1:n_phases
                        x = obj.FracturesNetwork.Fractures(f).State.Properties([key, num2str(i)]).Value;
                        x_Global(Start:End, i) = x;
                    end
                end
            end
        end
        function x_Global = CreateGlobalSinglePhaseVariables(obj, FineGrid, key)
            % For now only for fluid temperature, we use this function
            % instead of the function defined above which is "CreateGlobalVariables"
            x_Global   = zeros(sum([FineGrid(1:end).N]), 1);
            %% Reservoir
            x = obj.Reservoir.State.Properties(key).Value;
            x_Global(1:FineGrid(1).N) = x;
            %% Fractures
            if obj.FracturesNetwork.Active
                End = FineGrid(1).N;
                for f=1:obj.FracturesNetwork.NumOfFrac
                     Start = End+1;
                     End = Start + FineGrid(1+f).N - 1;
                     x = obj.FracturesNetwork.Fractures(f).State.Properties(key).Value;
                     x_Global(Start:End) = x;
                end
            end
        end
    end
end