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
        function ComputeFirstResidualNorm(obj, Residual, DiscretizationModel, LinearSolver)
            Nt_ADM = DiscretizationModel.ADMGrid.Ntot;
            % Get ADM Operators
            ADMRest = LinearSolver.R;
            % Restrict first residual
            Residual_ADM = ADMRest * Residual ./ sum(ADMRest, 2);
            obj.FirstResidual = Residual_ADM;
            % Compute Norms
            obj.FirstResidualNorm = zeros(obj.NumberOfEq,1);
            for eq = 1 : obj.NumberOfEq-1
                obj.FirstResidualNorm(eq) = norm(obj.FirstResidual((eq-1)*Nt_ADM+1:eq*Nt_ADM), 2);
            end
            obj.FirstResidualNorm(end) = norm(obj.FirstResidual(eq*Nt_ADM+1:end), 2);
            obj.NormCalculator.FirstResidualNorm = obj.FirstResidualNorm;
        end
        function converged = Check(obj, iter, residual, delta, Formulation, DiscretizationModel, State, LinearSolver)
            % Initialize
            converged = 0;
            Nt = DiscretizationModel.N;
            Nt_ADM = DiscretizationModel.ADMGrid.Ntot;
            
            ADMRest = LinearSolver.R;  % Get ADM Operators
            Residual_ADM = ADMRest * residual ./ sum(ADMRest, 2); % Restrict each residual and divide by number of cells in each coarse node
            
            % Compute Norms
            [ResidualNorm, equilibrium] =  obj.NormCalculator.CalculateResidualNorm(Residual_ADM, Nt_ADM, Formulation);
            [dp, dS] = obj.NormCalculator.CalculateSolutionNorm(delta, Nt, State);
            
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