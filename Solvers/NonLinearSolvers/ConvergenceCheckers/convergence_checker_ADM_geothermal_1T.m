% Convergence checker for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 2 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_ADM_geothermal_1T < convergence_checker_FS_geothermal_1T
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
            Nt = DiscretizationModel.N;
            Nt_ADM = DiscretizationModel.ADMGrid.Ntot;
            
            ADMRest = LinearSolver.R;  % Get ADM Operators
            Residual_ADM = ADMRest * residual ./ sum(ADMRest, 2); % Restrict each residual and divide by number of cells in each coarse node

            % Compute Norms            
            [ResidualNorm] = obj.NormCalculator.CalculateResidualNorm(Residual_ADM, Nt_ADM, Formulation);
            [dp, dT] = obj.NormCalculator.CalculateSolutionNorm(delta, Nt, State);
            
            converged = obj.CheckConvergenceCondition(iter,ResidualNorm,dp,dT);
        end
    end
end