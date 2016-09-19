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
        function SaveInitialState(obj, InitialState, Formulation)
            obj.State = copy(InitialState);
            Formulation.SavePhaseState();
        end
        function SetInitalGuess(obj, ProductionSystem)
            ProductionSystem.Reservoir.State = copy(obj.State);
        end
    end
end