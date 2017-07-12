%  FIM System Builder for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 14 October 2016
%Last modified: 14 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef fim_system_builder_ADM < fim_system_builder
    properties
 
    end
    methods
        function delta = UpdateState(obj, delta, ProductionSystem, Formulation, FluidModel, DiscretizationModel)
            % Update Reservoir State
            % Formulation.UpdatePandS(delta, ProductionSystem.Reservoir.State);
            
            % Solve FS equilibrium
            % Formulation.UpdateCompositions(ProductionSystem.Reservoir.State, FluidModel, DiscretizationModel.ReservoirGrid.N);
            % Update Reservoir State
            delta = Formulation.UpdateState(delta, ProductionSystem, FluidModel, DiscretizationModel);
            % UpdateWells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end
    end
end