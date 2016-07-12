% FIM coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 8 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FIM_Strategy < Coupling_Strategy
properties
    NLSolver
    TimeStepSelector
end
methods
    function obj = FIM_Strategy(name, NONLinearSolver)
        obj@Coupling_Strategy(name);
        obj.NLSolver = NONLinearSolver;
    end
    function [ProductionSystem, dt, Summary] = SolveTimeStep(obj, ProductionSystem, Formulation, maxDt, Summary)
        dt = objTimeStepSelector(maxDt);
        obj.NLSolver.Solve(ProductionSystem, Formulation, dt);
    end
end
end