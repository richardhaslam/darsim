% NL solver FS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 October 2016
%Last modified: 13 October 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef NL_Solver_FS < NL_Solver
    methods
        function UpdateState(obj, ProductionSystem, Formulation, FluidModel)
            obj.SystemBuilder.UpdateState(obj.Delta, ProductionSystem, Formulation, FluidModel);
        end
    end
end