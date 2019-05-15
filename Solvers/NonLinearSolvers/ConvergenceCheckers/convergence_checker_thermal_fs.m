% Convergence checker for Thermal simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_thermal_fs < convergence_checker
    properties
    end
    methods 
        function PrintTitles(obj, Residual)
            disp(['Initial residual norm: ', num2str(norm(Residual, inf))]);
            disp('');
            disp('        ||Residual||   ||delta p||   ||delta T||');
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            
            % Initialize
            converged = 0;
            
            
            disp(['Iter ' num2str(iter, '%02d') '    ' num2str(residual, '%5.5e'), '    ', num2str(dp, '%5.5e'), '    ', num2str(dS, '%5.5e')]);
            
            
            % Check convergence
            if (residula < obj.Tol) 
                converged = 1;
            end
        end
    end
end