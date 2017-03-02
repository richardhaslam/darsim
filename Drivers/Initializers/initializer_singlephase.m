% Initializer base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 8 November 2016
%Last modified: 8 November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef initializer_singlephase < initializer
    properties
        NLSolver
    end
    methods
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started Hydrostatic initialization');
            
            % Define initial values
            P_init = ones(DiscretizationModel.ReservoirGrid.N, 1)*1;
            
            % 1. Assign initial valus
            ProductionSystem.Reservoir.State.p = P_init;
            ProductionSystem.Reservoir.State.z(:,1) = 1;
            ProductionSystem.Reservoir.State.S(:,1) = 1;
            
            % 2. Update Composition of the phases (Flash)
            Formulation.SinglePhase = FluidModel.Flash(ProductionSystem.Reservoir.State);
            
            % 3 Compute Phase Density
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
        end
    end
end