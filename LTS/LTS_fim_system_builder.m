% LTS FIM System Builder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 12 April 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef LTS_fim_system_builder < fim_system_builder
    properties 
        LTSBCEnforcer = LTS_bc_enforcer_fim();
    end
    methods
        function Residual = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt)
            % Compute full residual (already with 0 transmissibilities where needed)
            Residual = Formulation.BuildResidual(ProductionSystem, DiscretizationModel, dt, obj.State);
            % Add b.c. to the residual
            Residual = obj.LTSBCEnforcer.AddBC2Residual(Residual, ProductionSystem, Formulation);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt)
            % Compute full Jacobian
            Jacobian = Formulation.BuildJacobian(ProductionSystem, DiscretizationModel, dt);
            % Add b.c. to the Jacobian
            Jacobian = obj.LTSBCEnforcer.AddBC2Jacobian(Jacobian, ProductionSystem, Formulation, DiscretizationModel);
        end
    end
end