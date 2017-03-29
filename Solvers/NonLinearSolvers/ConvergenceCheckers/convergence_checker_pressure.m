% Pressure convergence checker
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 26 July 2016
%Last modified: 26 July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_pressure < convergence_checker
    properties

    end
    methods 
        function PrintTitles(obj, Residual)
            disp(['Initial residual norm: ', num2str(norm(Residual, inf))]);
            disp('');
            disp('        ||Residual||   ||delta p||');
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State)
            % Initialize
            converged = 0;
            % Compute Norms
            Norm1 =  norm(residual, inf);
            Norm2 = norm(delta, inf);
            
            disp(['Iter ' num2str(iter) '    ' num2str(Norm1, '%5.5e'), '    ', num2str(Norm2,'%5.5e')]);
            
            %Check convergence
            if (Norm1 < obj.Tol) 
                converged = 1;
            end
        end
    end
end