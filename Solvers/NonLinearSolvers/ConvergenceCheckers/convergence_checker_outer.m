% Convergence checker for outer loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_outer < convergence_checker
    properties
    end
    methods 
        function Converged = Check(obj, State, State_old)
            Converged = 0;
            Norm1 = norm((State.p - State_old.p)/max(State.p), inf);
            Norm2 = norm(State.S - State_old.S, inf);
            
            if (Norm1 < obj.Tol && Norm2 < obj.Tol)
                Converged = 1;
            end
        end
    end
end