% Convergence checker for ADM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DARSim 2 Reservoir Simulator
%Author: Matteo Cusini
%TU Delft
%Created: 13 July 2016
%Last modified: 2 August 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef convergence_checker_ADM_geothermal_singlephase < convergence_checker_finescale_geothermal_singlephase
    properties
        OperatorsAssembler
    end
    methods
        function ComputeFirstResidualNorm(obj, Residual, RHS, DiscretizationModel, LinearSolver)
            Nt_ADM = DiscretizationModel.ADMGrid.Ntot;
            % Get ADM Operators
            ADMRest = LinearSolver.R;
            % Restrict first residual and first RHS
            obj.FirstResidual = ADMRest * Residual ./ sum(ADMRest, 2);
            obj.FirstRHS      = ADMRest * RHS      ./ sum(ADMRest, 2);
            
            % Compute Norms
            obj.FirstResidualNorm = zeros(obj.NumberOfEq,1);
            for eq = 1 : obj.NumberOfEq
                obj.FirstResidualNorm(eq) = norm(obj.FirstResidual((eq-1)*Nt_ADM+1:eq*Nt_ADM), 2);
                obj.FirstRHSNorm(eq)      = norm(obj.FirstRHS(     (eq-1)*Nt_ADM+1:eq*Nt_ADM), 2);
            end
            obj.NormCalculator.FirstResidualNorm = obj.FirstResidualNorm;
            obj.NormCalculator.FirstRHSNorm      = obj.FirstRHSNorm;
            
            % Initializing the "ResidualNorm" as empty, to make it clean for the coming iteration loop.
            obj.ResidualNorm = [];
            obj.RHSNorm = [];
        end
        function converged = Check(obj, iter, residual, RHS, delta, Formulation, DiscretizationModel, State, LinearSolver)
            Nt = DiscretizationModel.N;
            Nt_ADM = DiscretizationModel.ADMGrid.Ntot;
            
            % Get ADM Operators
            ADMRest = LinearSolver.R;
            Residual_ADM = ADMRest * residual ./ sum(ADMRest, 2); % Restrict each residual and divide by number of cells in each coarse node
            RHS_ADM      = ADMRest * RHS      ./ sum(ADMRest, 2); % Restrict each residual and divide by number of cells in each coarse node
            
            % Compute Norms
            [ obj.ResidualNorm(iter,:), obj.RHSNorm(iter,:) ] =  obj.NormCalculator.CalculateResidualNorm(Residual_ADM, RHS_ADM, Nt_ADM, Formulation);
            [dp, dT] = obj.NormCalculator.CalculateSolutionNorm(delta, Nt, State);
            
            obj.ResidualNorm( imag(obj.ResidualNorm) ~= 0 ) = NaN;
            obj.RHSNorm     ( imag(obj.RHSNorm     ) ~= 0 ) = NaN;
            dp              ( imag(dp              ) ~= 0 ) = NaN;
            dT              ( imag(dT              ) ~= 0 ) = NaN;
            
            converged = obj.CheckConvergenceCondition(iter,dp,dT);
        end
    end
end