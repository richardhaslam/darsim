% Coupling strategy base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Coupling_Strategy < handle
    properties
        Name
        TimeStepSelector
        Converged
    end
    methods
        function obj = Coupling_Strategy(name)
            obj.Name = name;
        end
    end
    methods (Abstract)
        [ProductionSystem, dt] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, maxDt)
    end
end