% Convergence checker for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 6 September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_ADM < convergence_checker_FS
    properties
    end
    methods
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State)
            
            % Initialize
            converged = 0;
            
            % Restrict Residual
            [R, ~] = DiscretizationModel.OperatorsHandler.AssembleFullOperators();
            residual_c = R * residual;
            
            % Compute Norms
            [massbalance, equilibrium] =  obj.NormCalculator.ResidualNorm(residual_c, DiscretizationModel.ADMGrid.N, Formulation);
            [dp, dS] = obj.NormCalculator.SolutionNorm(delta, DiscretizationModel.ADMGrid.N, State);
            
            disp(['Iter ' num2str(iter) '    ' num2str(massbalance, '%5.5e'), '    ', num2str(equilibrium,'%5.5e'), '    ', num2str(dp, '%5.5e'), '    ', num2str(dS, '%5.5e')]);
            
            %Check convergence
            if (massbalance < obj.Tol && equilibrium < obj.Tol && dp < obj.Tol * 1e2 && dS < obj.Tol * 1e2)
                converged = 1;
            end
        end
    end
end