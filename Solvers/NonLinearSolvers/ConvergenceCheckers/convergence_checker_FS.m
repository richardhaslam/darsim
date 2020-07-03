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
        function PrintTitles(obj)
            disp(['Initial Eq. 1 Residual Norm: ', num2str(obj.FirstResidualNorm(1), '%5.5e')]);
            disp(['Initial Eq. 2 Residual Norm: ', num2str(obj.FirstResidualNorm(2), '%5.5e')]);
            disp('');
            disp('           ||Residual Eq1||    ||Residual Eq2||    ||delta p||    ||delta S||');
        end
        function converged = Check(obj, iter, residual, RHS, delta, Formulation, DiscretizationModel, State, LinearSolver)
            % Initialize
            converged = 0;
            Nt = DiscretizationModel.N;
            
            % Compute Norms
            [ obj.ResidualNorm(iter,:), obj.RHSNorm(iter,:) ] = obj.NormCalculator.CalculateResidualNorm(residual, RHS, Nt, Formulation);
            [dp, dS] = obj.NormCalculator.CalculateSolutionNorm(delta, DiscretizationModel.N, State);
            
            disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(obj.ResidualNorm(iter,1), '%5.5e'), '        ' ...
                                                         , num2str(obj.ResidualNorm(iter,2), '%5.5e'), '        ' ...
                                                         , num2str(dp, '%5.5e'), '    ', num2str(dS, '%5.5e')]);
            
            % check if is stagnating
            stagnating = obj.Stagnating(obj.ResidualNorm(iter,:)./obj.FirstResidualNorm);
            
            %Check convergence
            if ( (obj.ResidualNorm(iter,1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)) || (obj.ResidualNorm(iter,1)/obj.RHSNorm(iter,1) < obj.ResidualTol(1)) ) && ...
               ( (obj.ResidualNorm(iter,2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)) || (obj.ResidualNorm(iter,2)/obj.RHSNorm(iter,2) < obj.ResidualTol(2)) ) && ...
               ( dp < obj.SolutionTol(1) && dS < obj.SolutionTol(2) )
                converged = 1;
            elseif iter>5 && dp < obj.SolutionTol(1) && dS < obj.SolutionTol(2) && ...
                   abs(obj.ResidualNorm(end-1,1) - obj.ResidualNorm(end-2,1)) < 1e-5 * obj.ResidualNorm(end-2,1) && ...
                   abs(obj.ResidualNorm(end  ,1) - obj.ResidualNorm(end-1,1)) < 1e-5 * obj.ResidualNorm(end  ,1) && ...
                   abs(obj.ResidualNorm(end-1,2) - obj.ResidualNorm(end-2,2)) < 1e-5 * obj.ResidualNorm(end-2,2) && ...
                   abs(obj.ResidualNorm(end  ,2) - obj.ResidualNorm(end-1,2)) < 1e-5 * obj.ResidualNorm(end  ,2)
               converged = 1;
            elseif (isnan(obj.ResidualNorm(iter,1)) || isnan(obj.ResidualNorm(iter,2)) || stagnating || isnan(dp) || isnan(dS))
                converged = -1;
            end
        end
    end
end