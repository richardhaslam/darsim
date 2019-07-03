% LTS Transport System Builder 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Ludovica Delpopolo
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTStransport_system_builder < transport_system_builder
    properties
        LTSBCEnforcer = LTS_bc_enforcer_seq();
    end
    methods
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt)
           Residual = Formulation.BuildTransportResidual(ProductionSystem, DiscretizationModel, dt, obj.State);
           % Add b.c. to the residual
          Residual = obj.LTSBCEnforcer.AddBC2Residual(Residual, ProductionSystem, Formulation, DiscretizationModel,  obj.State, dt);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt)
            Jacobian = Formulation.BuildTransportJacobian(ProductionSystem, DiscretizationModel, dt);
            % Add b.c. to the Jacobian
            Jacobian = obj.LTSBCEnforcer.AddBC2Jacobian(Jacobian, ProductionSystem, Formulation, DiscretizationModel, dt);
        end
        function SynchronizeProperties(obj, ProductionSystem, State_global, CellsSelected)
           Names = obj.State.Properties.keys;
            N_prop = double(obj.State.Properties.Count);
            % First the Reservoir
            for i=1:N_prop
                Index.Start = 1;
                Index.End = size(ProductionSystem.Reservoir.K,1);
                temp = ProductionSystem.Reservoir.State.Properties(Names{i});
                temp.Value = temp.Value .* CellsSelected.ActCells + State_global.Properties(Names{i}).Value .* (1 - CellsSelected.ActCells);
                % Save fractures state
                % To be implemented
            end
        end
    end
end