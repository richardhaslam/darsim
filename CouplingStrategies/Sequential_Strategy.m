% Sequential coupling strategy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 4 July 2016
%Last modified: 23 July 2016
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
    function AddPressureSolver(obj, pressuresolver)
        obj.PressureSolver = pressuresolver;
    end
    function dt = SolveTimeStep(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, maxDt)
        % This is the outer loop
        obj.itCount = 0;
        obj.Converged = 0;
        % Pressure timestep
        dt = obj.TimeStepSelector.ChooseTimeStep();
        % Phase Mobilities and total Mobility
        Formulation.ComputeTotalMobility(ProductionSystem, FluidModel);
        while obj.Converged == 0 && obj.itCount < obj.MaxIter
            % Solve pressure equation
            disp('Pressure Solver');
            disp('...............................................');
            obj.PressureSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            Formulation.ComputeFluxes(ProductionSystem, DiscretizationModel);
            disp('...............................................');
            % Check that velocity field is conservative
            conservative = Formulation.CheckMassConservation(DiscretizationModel.ReservoirGrid);
            if ~conservative
                return
            end
            % Choose stable timestep
            dt = obj.TimeStepSelector.StableTimeStep(ProductionSystem, DiscretizationModel, FluidModel, Formulation.U);
            % Solve transport
            disp('Transport Solver');
            disp('...............................................');
            obj.TransportSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
            disp('...............................................');
        end
        % Compute phase fluxes after converged solution
        ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
    end
end
end