% Convergence checker for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 6 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_ADM < convergence_checker
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
            
            % Restrict Residual
            [R, ~] = DiscretizationModel.OperatorsHandler.AssembleFullOperators();
            residual_c = R * residual;
            
            % Compute Norms
            Norm1 =  norm(residual_c, inf);
            Norm2 = norm(delta(1:N), inf)/max(State.p);
            Norm3 = norm(delta(N+1:end), inf);
            
            disp(['Iter ' num2str(iter) '    ' num2str(Norm1, '%5.5e'), '    ', num2str(Norm2,'%5.5e'), '    ', num2str(Norm3, '%5.5e')]);
            
            %Check convergence
            if (Norm1 < obj.Tol && Norm3 < obj.Tol)
                converged = 1;
            end
        end
    end
end