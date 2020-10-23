%  FIM System Builder base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 19 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%eh %
classdef fim_system_builder < system_builder
    properties
        NumberOfEq  
    end
    methods
        function ComputePropertiesAndDerivatives(obj, Formulation, ProductionSystem, FluidModel, DiscretizationModel)
            Formulation.ComputePropertiesAndDerivatives(ProductionSystem, FluidModel);
            Formulation.UpWindAndPhaseRockFluxes(DiscretizationModel, FluidModel.Phases, ProductionSystem);
        end
        function [Residual, RHS] = BuildResidual(obj, ProductionSystem, DiscretizationModel, Formulation, dt)
           [Residual, RHS] = Formulation.BuildFullResidual(ProductionSystem, DiscretizationModel, dt, obj.State);
        end
        function Jacobian = BuildJacobian(obj, ProductionSystem, Formulation, DiscretizationModel, dt)
            Jacobian = Formulation.BuildFullJacobian(ProductionSystem, DiscretizationModel, dt);
        end
        function SetUpSolutionChopper(obj, SolutionChopper, Formulation, ProductionSystem, DiscretizationModel)
            x = Formulation.GetPrimaryUnknowns(ProductionSystem, DiscretizationModel);
            SolutionChopper.DefineMaxDelta(x);
        end
        function delta = UpdateState(obj, delta, ProductionSystem, Formulation, FluidModel, DiscretizationModel)
            % Update State
            delta = Formulation.UpdateState(delta, ProductionSystem, FluidModel, DiscretizationModel);
            % Update Wells
            ProductionSystem.Wells.UpdateState(ProductionSystem.Reservoir, FluidModel);
        end
    end
end