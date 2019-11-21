% Convergence checker for Thermal simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Rhadityo
%TU Delft
%Created: 24 January 2018
%Last modified: 24 January 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_FS_geothermal_1T < convergence_checker
    properties
    end
    methods
        function ComputeFirstResidualNorm(obj, Residual, DiscretizationModel, LinearSolver)
            Nt = DiscretizationModel.N; 
            obj.FirstResidual = Residual;
            % Compute Norms
            obj.FirstResidualNorm = zeros(obj.NumberOfEq,1);
            for eq = 1 : obj.NumberOfEq-1
                obj.FirstResidualNorm(eq) = norm(obj.FirstResidual((eq-1)*Nt+1:eq*Nt), 2);
            end
            obj.FirstResidualNorm(end) = norm(obj.FirstResidual(eq*Nt+1:end), 2);
            obj.NormCalculator.FirstResidualNorm = obj.FirstResidualNorm;
        end
        function PrintTitles(obj)
            disp(['Initial Pressure    Residual Norm: ', num2str(obj.FirstResidualNorm(1), '%5.5e')]);
            disp(['Initial Temperature Residual Norm: ', num2str(obj.FirstResidualNorm(2), '%5.5e')]);
            disp('');
            disp('           ||Residual P||   ||Residual T||     ||delta P||    ||delta T||');
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            Nt = DiscretizationModel.N;

            % Compute Norms
            [ResidualNorm] = obj.NormCalculator.CalculateResidualNorm(residual, Nt, Formulation);
            [dp, dT] = obj.NormCalculator.CalculateSolutionNorm(delta, DiscretizationModel.N, State);
            
            converged = obj.CheckConvergenceCondition(iter,ResidualNorm,dp,dT);

        end
        function converged = CheckConvergenceCondition(obj,iter,ResidualNorm,dp,dT)
            disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(ResidualNorm(1), '%5.5e'), '      ' ...
                                                         , num2str(ResidualNorm(2), '%5.5e'), '      ' ...
                                                         , num2str(dp, '%5.5e'), '    ', num2str(dT, '%5.5e')]);
            % Initialize
            converged = 0;
            
            % check if it is stagnating
            stagnating = obj.Stagnating(ResidualNorm./obj.FirstResidualNorm);
            
            % Check convergence
            if ( ResidualNorm(1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1) && dp < obj.SolutionTol(1) ) && ...
               ( ResidualNorm(2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2) && dT < obj.SolutionTol(2) )
                converged = 1;
            elseif (dp < obj.SolutionTol(1)*1e-2) && (dT < obj.SolutionTol(2)*1e-2)
                converged = 1;
            elseif( isnan(ResidualNorm(1)) || isnan(ResidualNorm(2)) || stagnating || isnan(dp) || isnan(dT) || ~isreal([dp;dT]) )
                converged = -1;
            end
        end
    end
end