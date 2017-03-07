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
        function SaveInitialState(obj, InitialState, Formulation)            
            obj.State.CopyProperties(InitialState);
            Formulation.SavePhaseState();
        end
        function SetInitalGuess(obj, ProductionSystem)
            ProductionSystem.Reservoir.State.CopyProperties(obj.State);
        end
    end
end