% Sequential coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 18 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Sequential_Strategy < Coupling_Strategy
properties
    PressureSolver
    TransportSolver
    TimeStepSelector
end
methods
    function obj = Sequential_Strategy(name)
        obj@Coupling_Strategy(name);
    end
    function [ProductionSystem, dt] = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, maxDt)
        ProductionSystem = obj.PressureSolver.solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation);
        obj.ComputeVelocityField();
        obj.CheckVelocity();
        dt = obj.TimeStepSelector.ChooseTimestep(maxDt);
        Status = obj.TransportSolver.solve(Status);
    end
end
end