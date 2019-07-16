% ADM Preconditioner 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef preconditioner < handle
    properties
        Solver
    end
    methods
        function SetUp(obj, ProductionSystem, DiscretizationModel)
          
        end
        function x = Solve(A, rhs)
            % Matlab simple direct solver
            x = A\rhs;
        end
    end
end