% Convergence checker for MultiPhase Geothermal simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Arjan Marelis 
%TU Delft
%Created: 15 April 2020
%Last modified: 15 April 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_finescale_geothermal_multiphase < convergence_checker
    properties
    end
    methods
        function PrintTitles(obj)
            disp(['Initial Mass   Balance Residual Norm: ', num2str(obj.FirstResidualNorm(1), '%5.5e')]);
            disp(['Initial Energy Balance Residual Norm: ', num2str(obj.FirstResidualNorm(2), '%5.5e')]);
            disp('');
            disp('           ||Residual MB||   ||Residual EB||     ||delta P||    ||delta H||');
        end
        function converged = Check(obj, iter, residual, RHS, delta, Formulation, DiscretizationModel, State, LinearSolver)
            Nt = DiscretizationModel.N;

            % Compute Norms
            [ obj.ResidualNorm(iter,:), obj.RHSNorm(iter,:) ] = obj.NormCalculator.CalculateResidualNorm(residual, RHS, Nt, Formulation);
            [dp, dh] = obj.NormCalculator.CalculateSolutionNorm(delta, DiscretizationModel.N, State);
            
            obj.ResidualNorm( imag(obj.ResidualNorm) ~= 0 ) = NaN;
            obj.RHSNorm     ( imag(obj.RHSNorm     ) ~= 0 ) = NaN;
            dp              ( imag(dp              ) ~= 0 ) = NaN;
            dh              ( imag(dh              ) ~= 0 ) = NaN;
            
            converged = obj.CheckConvergenceCondition(iter,dp,dh);
        end
        function converged = CheckConvergenceCondition(obj,iter,dp,dh)
            disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(obj.ResidualNorm(iter,1), '%5.5e'), '      ' ...
                                                         , num2str(obj.ResidualNorm(iter,2), '%5.5e'), '      ' ...
                                                         , num2str(dp, '%5.5e'), '    ', num2str(dh, '%5.5e')]);
            % Initialize
            converged = 0;
            
            % check if it is stagnating
            stagnating = obj.Stagnating(obj.ResidualNorm(iter,:)./obj.FirstResidualNorm);
            
            % Check convergence
            if ( (obj.ResidualNorm(iter,1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.RHSNorm(iter,1) < obj.ResidualTol(1)) ) && ...
               ( (obj.ResidualNorm(iter,2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.RHSNorm(iter,2) < obj.ResidualTol(2)) ) && ...
               ( dp < obj.SolutionTol(1) && dh < obj.SolutionTol(2) )
                converged = 1;
            elseif iter>5 && dp < obj.SolutionTol(1) && dh < obj.SolutionTol(2) && ...
                   abs(obj.ResidualNorm(end-1,1) - obj.ResidualNorm(end-2,1)) < 1e-3 * obj.ResidualNorm(end-2,1) && ...
                   abs(obj.ResidualNorm(end  ,1) - obj.ResidualNorm(end-1,1)) < 1e-3 * obj.ResidualNorm(end  ,1) && ...
                   abs(obj.ResidualNorm(end-1,2) - obj.ResidualNorm(end-2,2)) < 1e-3 * obj.ResidualNorm(end-2,2) && ...
                   abs(obj.ResidualNorm(end  ,2) - obj.ResidualNorm(end-1,2)) < 1e-3 * obj.ResidualNorm(end  ,2)
               converged = 1;
            elseif (isnan(obj.ResidualNorm(iter,1)) || isnan(obj.ResidualNorm(iter,2)) || stagnating || isnan(dp) || isnan(dh))
                converged = -1;
            end
        end
    end
end