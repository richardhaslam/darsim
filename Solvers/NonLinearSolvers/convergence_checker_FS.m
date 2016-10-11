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
            converged = 0;
            % Compute Norms
            [massbalance, equilibrium] =  NormCalculator.ResidualNorm(residual, DiscretizationModel, Formulation);
            [dp, ds] = NormCalculator.SolutionNorm(delta);
            
            disp(['Iter ' num2str(iter) '    ' num2str(massbalance, '%5.5e'), '    ', num2str(equilibrium,'%5.5e'), '    ', num2str(dp, '%5.5e'), '    ', num2str(ds, '%5.5e')]);
            
            %Check convergence
            if (massbalance < obj.Tol && equilibrium < obj.Tol && dp < obj.Tol && ds < obj.Tol)
                converged = 1;
            end
        end
    end
end