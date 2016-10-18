% NL solver ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 October 2016
%Last modified: 13 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NL_Solver_ADM < NL_Solver
    methods
        function UpdateState(obj, ProductionSystem, Formulation, FluidModel)
            % Update P and S
            Formulation.UpdateState_ADM(obj.Delta, ProductionSystem, Formulation, FluidModel);
            Formulation.UpdateState
        end
    end
end