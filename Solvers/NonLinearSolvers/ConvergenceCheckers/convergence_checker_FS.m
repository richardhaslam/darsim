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
        FirstResidualNorm
        FirstResidual
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
            disp(['Initial Pressure   Residual Norm: ', num2str(obj.FirstResidualNorm(1), '%5.5e')]);
            disp(['Initial Saturation Residual Norm: ', num2str(obj.FirstResidualNorm(2), '%5.5e')]);
            disp('');
            disp('           ||Residual P||    ||Residual S||    ||delta p||    ||delta S||');
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            % Initialize
            converged = 0;
            Nt = DiscretizationModel.N;
            
            % Compute Norms
            [ResidualNorm, equilibrium] = obj.NormCalculator.CalculateResidualNorm(residual, Nt, Formulation);
            [dp, dS] = obj.NormCalculator.CalculateSolutionNorm(delta, DiscretizationModel.N, State);
            
            disp(['Iter ', num2str(iter, '%02d') '-->   ', num2str(ResidualNorm(1)/obj.FirstResidualNorm(1), '%5.5e'), '      ' ...
                                                         , num2str(ResidualNorm(2)/obj.FirstResidualNorm(2), '%5.5e'), '      ' ...
                                                         , num2str(dp, '%5.5e'), '    ', num2str(dS, '%5.5e')]);
            
            % check if is stagnating
            stagnating = obj.Stagnating(ResidualNorm./obj.FirstResidualNorm);
            
            %Check convergence
            if ( (ResidualNorm(1) < obj.ResidualTol(1)) || (ResidualNorm(1)/obj.FirstResidualNorm(1) < obj.ResidualTol(1)) ) && ...
               ( (ResidualNorm(2) < obj.ResidualTol(2)) || (ResidualNorm(2)/obj.FirstResidualNorm(2) < obj.ResidualTol(2)) ) && ...
               ( dp < obj.SolutionTol(1) && dS < obj.SolutionTol(2) )
                converged = 1;
            elseif (isnan(ResidualNorm(1)) || isnan(ResidualNorm(2)) || stagnating || isnan(dp) || isnan(dS))
                converged = -1;
            end
        end
    end
end