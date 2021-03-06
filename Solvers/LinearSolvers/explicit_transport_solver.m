%  Explicit transport solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef explicit_transport_solver < handle
    properties
        itCount
        SystemBuilder
        Converged;
    end
    methods
        function obj = explicit_transport_solver()
            obj.itCount = 1;
        end
        function SetUp(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel, dt)
            % Virtual call
        end
        function SetUpLinearSolver(obj, ProductionSystem, DiscretizationModel)
            % Virtual call
        end
        function Solve(obj, ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt)
            disp('explicit transport solver');      
            % 1. Compute fractional flow
            obj.SystemBuilder.ComputePropertiesAndDerivatives(Formulation, ProductionSystem, FluidModel, DiscretizationModel);
            
            % 2. Assemble system matrix
            Formulation.ViscousMatrix(DiscretizationModel.ReservoirGrid);
            
            % 3. Solve explicit transport
            Formulation.UpdateSaturationExplicitly(ProductionSystem, DiscretizationModel, dt);
            
            obj.Converged = 1;
        end
    end
end