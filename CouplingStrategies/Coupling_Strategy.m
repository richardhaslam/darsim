% Coupling strategy base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 8 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Coupling_Strategy
properties
    Name
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