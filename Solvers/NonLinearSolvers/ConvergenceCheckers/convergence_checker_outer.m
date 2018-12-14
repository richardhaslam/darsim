% Convergence checker for outer loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_outer < convergence_checker
    properties
    end
    methods 
        function Converged = Check(obj, State, State_old)
            Converged = 0;
            p = State.Properties('P_1').Value;
            p_old = State_old.Properties('P_1').Value;
            s = State.Properties('S_1').Value;
            s_old = State_old.Properties('S_1').Value;
            Norm1 = norm((p - p_old)/max(p), inf);
            Norm2 = norm(s - s_old, inf);
            
            if (Norm1 < obj.SolutionTol && Norm2 < obj.SolutionTol)
                Converged = 1;
            end
        end
    end
end