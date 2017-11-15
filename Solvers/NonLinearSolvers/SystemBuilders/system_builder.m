%  System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef system_builder < handle
    properties
        State
    end
    methods (Abstract)
        obj = ComputePropertiesAndDerivatives(obj);
        obj = BuildResidual(obj);
        obj = BuildJacobian(obj);
    end
    methods
        function obj = system_builder()
            obj.State = status();
        end
        function SaveInitialState(obj, ProductionSystem, Formulation)
            % First the Reservoir
            obj.State.CopyProperties(ProductionSystem.Reservoir.State);
            Formulation.SavePhaseState();
            
            % Save fractures state
            for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                Names = obj.State.Properties.keys;
                N_prop = double(obj.State.Properties.Count);
                for i = 1:N_prop
                    temp = obj.State.Properties(Names{i});
                    temp.Value = [obj.State.Properties(Names{i}).Value; ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(Names{i}).Value];
                end
            end
        end
        function SetInitalGuess(obj, ProductionSystem)
            Names = obj.State.Properties.keys;
            N_prop = double(obj.State.Properties.Count);
            % First the Reservoir
            for i=1:N_prop
                Index.Start = 1;
                Index.End = size(ProductionSystem.Reservoir.K,1);
                temp = ProductionSystem.Reservoir.State.Properties(Names{i});
                temp.Value = obj.State.Properties(Names{i}).Value(Index.Start:Index.End,:);
                % Save fractures state
                for f = 1:ProductionSystem.FracturesNetwork.NumOfFrac
                    Index.Start = Index.End+1;
                    Index.End = Index.Start + size(ProductionSystem.FracturesNetwork.Fractures(f).K,1)-1;
                    temp = ProductionSystem.FracturesNetwork.Fractures(f).State.Properties(Names{i});
                    temp.Value = obj.State.Properties(Names{i}).Value(Index.Start:Index.End,:);
                end
            end
        end
    end
end