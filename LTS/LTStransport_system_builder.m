% LTS Transport System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTStransport_system_builder < transport_system_builder
    properties
    end
    methods
        function ComputePropertiesAndDerivatives(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel)
            Formulation.UpdateFractionalFlow(ProductionSystem, FluidModel);
            Formulation.ComputeDfDS(ProductionSystem, FluidModel);
        end
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt, CellsSelected)
           Residual = Formulation.BuildTransportResidualLTS(ProductionSystem, DiscretizationModel, dt, obj.State, CellsSelected);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt, CellsSelected)
            Jacobian = Formulation.BuildTransportJacobianLTS(ProductionSystem, DiscretizationModel, dt, CellsSelected);
        end
        function SynchronizeProperties(obj, ProductionSystem, State_global, CellsSelected)
           Names = obj.State.Properties.keys;
            N_prop = double(obj.State.Properties.Count);
            % First the Reservoir
            for i=1:N_prop
                Index.Start = 1;
                Index.End = size(ProductionSystem.Reservoir.K,1);
                temp = ProductionSystem.Reservoir.State.Properties(Names{i});
                temp.Value = temp.Value .* CellsSelected.ActCells + ...
                             State_global.Properties(Names{i}).Value(Index.Start:Index.End,:) .* (1 - CellsSelected.ActCells);
                % Save fractures state
                % To be implement
            end
        end
    end
end