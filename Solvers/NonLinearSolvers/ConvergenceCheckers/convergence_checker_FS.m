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
            disp('        ||Residual||   ||Equilibrium||   ||delta p||   ||delta S||');
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            
            % Initialize
            converged = 0;
            % Compute Norms
            [massbalance, equilibrium] =  obj.NormCalculator.ResidualNorm(residual, DiscretizationModel.ReservoirGrid.N, Formulation);
            [dp, dS] = obj.NormCalculator.SolutionNorm(delta, DiscretizationModel.N, State);
            
            disp(['Iter ' num2str(iter, '%02d') '    ' num2str(massbalance, '%5.5e'), '    ', num2str(equilibrium,'%5.5e'), '    ', num2str(dp, '%5.5e'), '    ', num2str(dS, '%5.5e')]);
            
            % check if is stagnating
            stagnating = obj.Stagnating(massbalance);
            
            %Check convergence
            if (massbalance < obj.Tol && equilibrium < obj.Tol && dp < obj.Tol * 1e2 && dS < obj.Tol * 1e2)
                converged = 1;
            elseif (isnan(massbalance) || stagnating || isnan(dp) || isnan(dS))
                converged = -1;
            end
        end
    end
end