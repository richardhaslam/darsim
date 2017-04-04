% ADM Preconditioner 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 30 March 2017
%Last modified: 30 March 2017
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