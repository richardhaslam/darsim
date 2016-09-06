% Convergence checker for FS NLSolver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 13 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_FS < convergence_checker
    properties
    end
    methods 
        function PrintTitles(obj, Residual)
            disp(['Initial residual norm: ', num2str(norm(Residual, inf))]);
            disp('');
            disp('        ||Residual||   ||delta p||   ||delta S||');
        end
        function converged = Check(obj, iter, residual, delta, DiscretizationModel, State)
            
            % Initialize
            N = DiscretizationModel.ReservoirGrid.N;
            converged = 0;
            % Compute Norms
            Norm1 =  norm(residual, inf);
            Norm2 = norm(delta(1:N), inf)/max(State.p);
            Norm3 = norm(delta(N+1:end), inf);
            
            disp(['Iter ' num2str(iter) '    ' num2str(Norm1, '%5.5e'), '    ', num2str(Norm2,'%5.5e'), '    ', num2str(Norm3, '%5.5e')]);
            
            %Check convergence
            if (Norm1 < obj.Tol && Norm2 < obj.Tol && Norm3 < obj.Tol)
                converged = 1;
            end
        end
    end
end