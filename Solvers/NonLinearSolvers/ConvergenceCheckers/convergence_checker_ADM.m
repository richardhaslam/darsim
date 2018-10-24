% Convergence checker for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_ADM < convergence_checker_FS
    properties
        OperatorsAssembler
    end
    methods
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            
            % Initialize
            converged = 0;
            
            % Restrict Residual
            % Get ADM Operators
            R = LinearSolver.R;
            residual_c = R * residual;
            residual_c = residual_c ./ sum(R, 2); % divide by number of cells in each coarse node
            [N_adm, Nf] = size(R);
            N_adm = N_adm/2;
            Nf = Nf/2;
            
            % Compute Norms
            [massbalance, equilibrium] =  obj.NormCalculator.ResidualNorm(residual_c, N_adm, Formulation);
            [dp, dS] = obj.NormCalculator.SolutionNorm(delta, Nf, State);
            
            disp(['Iter ' num2str(iter) '    ' num2str(massbalance, '%5.5e'), '    ', num2str(equilibrium,'%5.5e'), '    ', num2str(dp, '%5.5e'), '    ', num2str(dS, '%5.5e')]);
            
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