% Incompressible pressure solver (Not used any more)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef incompressible_pressure_solver < handle
    properties
        A   % Pressure system matrix
        rhs % right-hand side
        LinearSolver % Linear Solver
    end
    methods
        function Solve(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt)
            disp('incompressible pressure solver');
            % 1. Compute Transmissibilities
            Formulation.ComputeTransmissibilities(ProductionSystem, DiscretizationModel);
            
            % 2. Construct matrix and rhs
            obj.A = Formulation.BuildIncompressiblePressureMatrix(DiscretizationModel);
            obj.rhs = Formulation.BuildIncompressibleRHS(ProductionSystem, DiscretizationModel, FluidModel);
            
            % 3. Add wells
            [obj.A, obj.rhs] = Formulation.AddWellsToPressureSystem(ProductionSystem.Wells, ProductionSystem.Reservoir.K, obj.A, obj.rhs);
            
            % 4. Linear Solver
            ProductionSystem.Reservoir.State.p = obj.LinearSolver.Solve(obj.A, obj.rhs);
        end
    end
end