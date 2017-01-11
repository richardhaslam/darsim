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
            P_init = ones(DiscretizationModel.ReservoirGrid.N, 1)*1.013e5;
            z_init = ones(DiscretizationModel.ReservoirGrid.N, 1)*0;
            
            
            % 1. Assign initial valus
            ProductionSystem.Reservoir.State.p = P_init;
            ProductionSystem.Reservoir.State.z(:,1) = z_init;
            ProductionSystem.Reservoir.State.z(:,2) = 1 - z_init;
           
            obj.Equilibrate(ProductionSystem, FluidModel, Formulation, DiscretizationModel);
        end
        function Equilibrate(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            equilibrium = 0;
            % Save initial State
            %obj.NLSolver.SystemBuilder.SaveInitialState(ProductionSystem.Reservoir.State, Formulation);
            P_init = ProductionSystem.Reservoir.State.p;
            while (~equilibrium)
                p0 = ProductionSystem.Reservoir.State.p;
                
                % 1. Update Composition of the phases (Flash)
                Formulation.SinglePhase = FluidModel.Flash(ProductionSystem.Reservoir.State);
                
                % 2 Compute Phase Density
                FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
                
                % 3. Update S based on components mass balance
                FluidModel.ComputePhaseSaturation(ProductionSystem.Reservoir.State, Formulation.SinglePhase);
                
                % 4. Total Density
                ProductionSystem.Reservoir.State.rhoT = FluidModel.ComputeTotalDensity(ProductionSystem.Reservoir.State.S, ProductionSystem.Reservoir.State.rho);
                
                % 5. Compute initial Pc
                ProductionSystem.Reservoir.State.pc = FluidModel.ComputePc(ProductionSystem.Reservoir.State.S);
                
                % 6. Compute pressure 
                ProductionSystem.Reservoir.State.p = P_init + ProductionSystem.Reservoir.State.rhoT .* Formulation.GravityModel.g .* (ProductionSystem.Reservoir.Depth - DiscretizationModel.ReservoirGrid.Depth);
                delta = ProductionSystem.Reservoir.State.p - p0;
                
                if (delta < 1e-3) % If solution doesn t change I am at equilibrium
                    equilibrium = 1;
                end
            end
        end
    end
end

