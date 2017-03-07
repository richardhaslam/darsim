% InitializerSinglePhase
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
        function obj = initializer_singlephase(names, values)
            obj@initializer(names, values);
        end
        function ComputeInitialState(obj, ProductionSystem, FluidModel, Formulation, DiscretizationModel)
            disp('Started single phase initialization');
            % 2. Update Composition of the phases (Flash)
            Formulation.SinglePhase = FluidModel.Flash(ProductionSystem.Reservoir.State);
            
            % 3 Compute Phase Density
            FluidModel.ComputePhaseDensities(ProductionSystem.Reservoir.State);
            
            % Output initial status:      
            disp('Initial conditions:')
            disp(['reservoir pressure:' num2str(max(ProductionSystem.Reservoir.State.Properties('P_1').Value/1e5)), ' bar']);
            disp(['Single phase reservoir']);
            disp(['reservoir temperature: ', num2str(ProductionSystem.Reservoir.Temp)]);
            disp('---------------------------------------------------------');
            disp(char(5));
        end
    end
end