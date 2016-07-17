% Linear solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 15 July 2016
%Last modified: 15 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef linear_solver < handle
    properties
        CondNumber
    end
    methods
        function x = Solve(obj, A, rhs)
            % Use Matlab direct solver
            obj.CondNumber = cond(A);
            x = A\rhs;
        end
    end
end