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
    MaxIter
    itCount
    Converged
end
methods
    function obj = Sequential_Strategy(name)
        obj@Coupling_Strategy(name);
    end
    function dt = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, maxDt)
        % This is the outer loop
        obj.itCount = 0;
        obj.Converged = 0;
        while obj.Converged == 0 && obj.itCount < obj.MaxIter
            % Solve pressure equation
            obj.PressureSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            Formulation.ComputeFluxes(ProductionSystem, DiscretizationModel);
            % Check that velocity field is conservative
            conservative = Formulation.CheckMassConservation();
            if ~conservative
                return
            end
            % Choose stable timestep
            dt = obj.TimeStepSelector.ChooseTimestep(maxDt);
            % Solve transport
            obj.TransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
        end
    end
end
end