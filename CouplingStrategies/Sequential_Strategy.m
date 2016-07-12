% Sequential coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 12 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Sequential_Strategy < Coupling_Strategy
properties
    PressureSolver
    TransportSolver
end
methods
    function obj = Sequential_Strategy(name)
        obj@Coupling_Strategy(name);
    end
    function Status = SolveTimeStep(obj, Status)
        Status = obj.PressureSolver.solve(Status);
        obj.ComputeVelocityField;
        obj.CheckVelocity();
        Status = obj.TransportSolver.solve(Status);
    end
end
end