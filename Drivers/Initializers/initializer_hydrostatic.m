% Initializer base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 8 November 2016
%Last modified: 8 November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer_hydrostatic < initializer
    properties
        NLSolver
    end
    methods
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started Hydrostatic initialization');
            
             % Define initial values
            P_init = ones(DiscretizationModel.ReservoirGrid.N, 1)*1e5;
            z_init = ones(DiscretizationModel.ReservoirGrid.N, 1)*0.05;
            
            
            % 1. Assign initial valus
            ProductionSystem.Reservoir.State.p = P_init;
            ProductionSystem.Reservoir.State.z(:,1) = z_init;
            ProductionSystem.Reservoir.State.z(:,2) = 1 - z_init;
           
            Formulation.SinglePhase = FluidModel.InitializeReservoir(ProductionSystem.Reservoir.State);
            
            obj.Equilibrate(ProductionSystem, FluidModel, Formulation, DiscretizationModel);
        end
        function Equilibrate(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            equilibrium = 0;
            % Save initial State
            obj.NLSolver.SystemBuilder.SaveInitialState(ProductionSystem.Reservoir.State, Formulation);
            while (~equilibrium)
                % do nn-linear solve
                obj.NLSolver.Solve(ProductionSystem, FluidModel, DiscretizationModel, Formulation, dt);
                if (delta < 1e-6) % If solution doesn t change I am at equilibrium
                    equilibrium = 1;
                end
            end
        end
    end
end

